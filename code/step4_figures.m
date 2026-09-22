function step4_figures(cfg)
% STEP4_FIGURES  Build the computed manuscript figures as vector PDFs in
% results/figures, reading the result files written by steps 1-3.
%   Fig2  dataset characterisation
%   Fig3  model ablation with 95% confidence intervals
%   Fig4  mechanism of the target-set-size confound; panel (c) is reserved
%         for the example subnetwork, which is drawn separately
%   FigS1 drug target-profile embedding
%   FigS2 feed-forward-loop convergence on shared effectors
% Fig5 (baseline) and Fig6 (external validation) are produced by their own
% steps. All plotting helpers are local functions below.

fprintf('=== step 4: figures ===\n');
BLUE=[0.165 0.471 0.839]; ORANGE=[0.922 0.408 0.204]; MUTED=[0.537 0.529 0.506];
GRAY=[0.722 0.718 0.690]; FILL=[0.917 0.945 0.988]; GREEN=[0.047 0.639 0.047]; RED=[0.890 0.286 0.282];
outdir = cfg.figures;

% -------------------- Figure 2: dataset characterisation -------------
G = load(fullfile(cfg.results,'graph.mat')); zip=G.zip; pairOK=G.pairOK; cl=G.cl;
f=figure('Position',[80 80 1200 360],'Color','w');
z=zip(pairOK & ~isnan(zip));
axa=subplot(1,3,1); hold on; lo=pctl(z,0.5); hi=pctl(z,99.5);
[fd,xi]=kde1(z,lo,hi,400); f0=interp1(xi,fd,0); f10=interp1(xi,fd,10);
mA=xi<0; mN=xi>=0 & xi<=10; mS=xi>10;
xa=[xi(mA) 0]; ya=[fd(mA) f0]; xn=[0 xi(mN) 10]; yn=[f0 fd(mN) f10]; xs=[10 xi(mS)]; ys=[f10 fd(mS)];
patch([xa fliplr(xa)],[ya zeros(size(ya))],BLUE,'EdgeColor','none','FaceAlpha',0.85);
patch([xn fliplr(xn)],[yn zeros(size(yn))],GRAY,'EdgeColor','none','FaceAlpha',0.85);
patch([xs fliplr(xs)],[ys zeros(size(ys))],ORANGE,'EdgeColor','none','FaceAlpha',0.9);
plot(xi,fd,'-','Color',[0.05 0.05 0.05],'LineWidth',1.4); yl=[0 max(fd)*1.22];
plot([0 0],yl,'--','Color',MUTED,'LineWidth',1); plot([10 10],yl,'-','Color',ORANGE,'LineWidth',1.4);
text(lo+0.04*(hi-lo),yl(2)*0.93,'antag.','Color',BLUE,'FontWeight','bold','FontSize',9);
text(4,yl(2)*0.5,'neutral','Color',[0.45 0.45 0.43],'FontWeight','bold','FontSize',9,'HorizontalAlignment','center');
text(12,yl(2)*0.85,'syn.','Color',ORANGE,'FontWeight','bold','FontSize',9);
xlim([lo hi]); ylim(yl); xlabel('ZIP synergy score'); ylabel('density'); title('(a) Synergy is rare'); fig_style(axa);
axb=subplot(1,3,2); hold on; axis equal off;
nSyn=sum(z>=10); nNon=sum(z<=0); nAmb=sum(z>0 & z<10); tot=nSyn+nNon+nAmb;
vals=[nSyn nAmb nNon]; dcols={ORANGE,GRAY,BLUE}; dlab={'synergistic','ambiguous','non-syn.'};
Ro=1; Ri=0.62; ang=pi/2;
for i=1:3
    a1=ang; a0=ang-(vals(i)/tot)*2*pi; th=linspace(a1,a0,90);
    patch([Ro*cos(th) fliplr(Ri*cos(th))],[Ro*sin(th) fliplr(Ri*sin(th))],dcols{i},'EdgeColor','w','LineWidth',1.2,'DisplayName',sprintf('%s (%.0f%%)',dlab{i},100*vals(i)/tot));
    am=(a1+a0)/2; if vals(i)/tot>0.05, text(0.8*cos(am),0.8*sin(am),sprintf('%.0f%%',100*vals(i)/tot),'HorizontalAlignment','center','FontSize',9); end
    ang=a0;
end
text(0,0,addcommas(tot),'HorizontalAlignment','center','FontSize',9,'Color',[0.32 0.32 0.30]);
title('(b) Synergy is imbalanced','FontWeight','bold','FontSize',13,'Color',[0.05 0.05 0.05]);
lg=legend('Location','southoutside'); set(lg,'Box','off','FontSize',8);
axc=subplot(1,3,3); hold on;
codes={'BR','CNS','CO','LE','ME','LC','OV','PR','RE'}; tname={'Breast','CNS','Colon','Leukemia','Melanoma','Lung','Ovarian','Prostate','Renal'};
cell2t=containers.Map('KeyType','char','ValueType','char'); exprfile='';
if isfield(cfg,'f_expr') && ~isempty(cfg.f_expr)
    if exist(cfg.f_expr,'file'), exprfile=cfg.f_expr;
    elseif isfield(cfg,'data'), exprfile=fullfile(cfg.data,cfg.f_expr); end
end
if ~isempty(exprfile) && exist(exprfile,'file')
    cn=read_cellminer(exprfile);
    for j=1:numel(cn), nmj=cn{j}; kk=strfind(nmj,':'); if isempty(kk), continue; end; code=nmj(1:kk(1)-1); idx=find(strcmp(codes,code),1); if isempty(idx), continue; end; cell2t(normcell(nmj))=tname{idx}; end
end
tc=containers.Map(tname,num2cell(zeros(1,numel(tname))));
for i=1:numel(cl), if ~pairOK(i)||isnan(zip(i))||zip(i)<10, continue; end; key=normcell(cl{i}); if isKey(cell2t,key), tt=cell2t(key); tc(tt)=tc(tt)+1; end; end
tv=cellfun(@(t)tc(t),tname); [tv,si]=sort(tv,'descend'); tl=tname(si); warm=bluemap(numel(tv));
for i=1:numel(tv), barh(numel(tv)-i+1,tv(i),0.72,'FaceColor',warm(i,:),'EdgeColor','none'); text(tv(i)+max(max(tv),1)*0.02,numel(tv)-i+1,num2str(tv(i)),'FontSize',8,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle'); end
set(axc,'YTick',1:numel(tv),'YTickLabel',fliplr(tl),'FontSize',8); ylim([0.4 numel(tv)+0.6]); if max(tv)>0, xlim([0 max(tv)*1.2]); end
xlabel('synergistic pairs'); title('(c) Synergy varies by tissue'); fig_style(axc); set(axc,'FontSize',8);
savepdf(f,fullfile(outdir,'Fig2_data'));

% -------------------- Figure 3: ablation with 95% CIs ----------------
M = load(fullfile(cfg.results,'model_results.mat'));
f=figure('Position',[80 80 1180 340],'Color','w');
xpos=1:4; off=0.16;
ax1=subplot(1,3,1); hold on; ablation_panel(ax1, M.om_roc, 0.5, 0.5, 1.0, '(a) AUROC by feature set','AUROC (95% CI)', BLUE, ORANGE);
ax2=subplot(1,3,2); hold on; ablation_panel(ax2, M.om_prc, M.om_base, 0, 0.55, '(b) AUPRC by feature set','AUPRC (95% CI)', BLUE, ORANGE);
ax3=subplot(1,3,3); hold on;
gi=M.gimp; v=[gi(3) gi(2) gi(1)]; bm=bluemap(3); [~,rk]=sort(v,'descend'); bc=zeros(3,3); bc(rk,:)=bm;
for i=1:3, barh(i,v(i),0.6,'FaceColor',bc(i,:),'EdgeColor','none'); text(v(i)+0.015,i,sprintf('%.2f',v(i)),'FontSize',10,'Color',[0.25 0.25 0.24]); end
set(ax3,'YTick',1:3,'YTickLabel',{'network','cell','drug'}); ylim([0.4 3.6]); xlim([0 max(v)+0.16]);
xlabel('importance share'); title('(c) Feature importance'); fig_style(ax3);
savepdf(f,fullfile(outdir,'Fig3_results'));

% -------------------- Figure 4: mechanism ----------------------------
% (a) apparent signal, (b) the target-set-size confound, (c) reserved for
% the example feed-forward-loop subnetwork drawn separately.
F = load(fullfile(cfg.results,'pair_features.mat'));
fn=F.featnames; yv=F.y; labm=~isnan(yv);
f=figure('Position',[80 80 1200 340],'Color','w');
axm1=subplot(1,3,1);
boxpts(axm1, F.X(labm,strcmp(fn,'convFFL')), yv(labm), BLUE, ORANGE, '(a) FFLs per pair','feed-forward loops');
axm2=subplot(1,3,2); hold(axm2,'on');
tcount = F.X(:,strcmp(fn,'nTA')) + F.X(:,strcmp(fn,'nTB'));
scatterconf(axm2, tcount(labm), F.X(labm,strcmp(fn,'convFFL')), yv(labm), BLUE, ORANGE);
title(axm2,'(b) FFLs track target-set size','FontWeight','bold','FontSize',13,'Color',[0.05 0.05 0.05]);
xlabel(axm2,'targets per pair'); ylabel(axm2,'feed-forward loops'); fig_style(axm2);
axm3=subplot(1,3,3); axis(axm3,'off');
title(axm3,'(c) Example FFL subnetwork','FontWeight','bold','FontSize',13,'Color',[0.05 0.05 0.05]);
savepdf(f,fullfile(outdir,'Fig4_mechanism'));

% -------------------- Supplementary: feature space and convergence ----
biology_panels(cfg, outdir, BLUE, ORANGE);
fprintf('wrote Fig2, Fig3, Fig4 and supplementary panels to %s\n', outdir);
end

% ===================== supplementary biology =========================
function biology_panels(cfg, outdir, BLUE, ORANGE)
% Two descriptive panels: the drug target-profile embedding, whose leading
% axis tracks target-set size, and the per-effector convergence of
% feed-forward loops in synergistic versus non-synergistic pairs.
    S = load(fullfile(cfg.results,'signor.mat'));
    G = load(fullfile(cfg.results,'graph.mat'));
    P = load(fullfile(cfg.results,'pair_features.mat'));
    dsym=S.dsym; Ddir=logical(S.Ddir);
    d2t = containers.Map(G.dkeys,G.dvals);
    sig = containers.Map('KeyType','char','ValueType','double');
    for i=1:numel(dsym), sig(upper(dsym{i}))=i; end

    % drug target-profile embedding
    drugs = unique([P.pairA(:); P.pairB(:)]); keep=false(numel(drugs),1);
    for i=1:numel(drugs), keep(i)=isKey(d2t,drugs{i}) && ~isempty(d2t(drugs{i})); end
    drugs=drugs(keep);
    allg={}; for i=1:numel(drugs), allg=[allg, d2t(drugs{i})]; end
    allg=unique(allg); gi=containers.Map(allg,num2cell(1:numel(allg)));
    DM=zeros(numel(drugs),numel(allg));
    for i=1:numel(drugs), g=d2t(drugs{i}); for k=1:numel(g), DM(i,gi(g{k}))=1; end, end
    ntar=sum(DM,2); Mc=DM-mean(DM,1); [U,Sg,~]=svd(Mc,'econ'); sc=U*Sg;
    ev=diag(Sg).^2; ev=ev/sum(ev);
    f=figure('Position',[80 80 520 460],'Color','w'); ax=axes('Parent',f); hold(ax,'on');
    scatter(ax, sc(:,1), sc(:,2), 46, log10(ntar), 'filled','MarkerEdgeColor','w','LineWidth',0.5);
    colormap(ax, viridis_local()); cb=colorbar(ax); ylabel(cb,'log_{10} targets');
    xlabel(ax,sprintf('PC1 (%.1f%%)',100*ev(1))); ylabel(ax,sprintf('PC2 (%.1f%%)',100*ev(2)));
    title(ax,sprintf('Drug target-profile space (n=%d)',numel(drugs)),'FontWeight','bold');
    fig_style(ax); savepdf(f,fullfile(outdir,'FigS1_drug_space'));

    % feed-forward-loop convergence on shared downstream effectors
    tg=@(k) tg_in_signor(d2t,sig,k);
    hs=reach_profile(P.pairA(P.y==1), P.pairB(P.y==1), tg, Ddir, numel(dsym));
    hn=reach_profile(P.pairA(P.y==0), P.pairB(P.y==0), tg, Ddir, numel(dsym));
    keepn=find((hs+hn)>0.15); dev=abs(hs-hn); lim=1.08*max([hs(keepn);hn(keepn)]);
    f=figure('Position',[80 80 520 500],'Color','w'); ax=axes('Parent',f); hold(ax,'on');
    plot(ax,[0 lim],[0 lim],'--','Color',[0.55 0.55 0.55],'LineWidth',1);
    for j=keepn(:)'
        if dev(j)>0.06, c=ORANGE; if hs(j)<hn(j), c=BLUE; end; else, c=[0.72 0.75 0.79]; end
        scatter(ax,hn(j),hs(j),60,c,'filled','MarkerEdgeColor','w','LineWidth',0.5);
        if dev(j)>0.08 || (hs(j)+hn(j))>0.9, text(ax,hn(j)+0.01*lim,hs(j)+0.01*lim,dsym{j},'FontSize',8.5); end
    end
    axis(ax,[0 lim 0 lim]); axis(ax,'square');
    xlabel(ax,'FFL reach per non-synergistic pair'); ylabel(ax,'FFL reach per synergistic pair');
    title(ax,'FFL convergence on effectors','FontWeight','bold'); fig_style(ax);
    savepdf(f,fullfile(outdir,'FigS2_convergence'));
end

function hc = reach_profile(A, B, tg, Ddir, N)
% Mean number of feed-forward loops per pair whose shared downstream node is
% each gene. Only present regulator edges are visited.
    hc=zeros(N,1); np=0;
    for p=1:numel(A)
        TA=tg(A{p}); TB=tg(B{p}); if isempty(TA)||isempty(TB), continue; end
        RA=Ddir(TA,:); RB=Ddir(TB,:);
        [ia,jb]=find(Ddir(TA,TB));
        for e=1:numel(ia)
            x=TA(ia(e)); yv=TB(jb(e)); if x==yv, continue; end
            z=find(RA(ia(e),:) & RB(jb(e),:)); z(z==x|z==yv)=[]; hc(z)=hc(z)+1;
        end
        [ib,ja]=find(Ddir(TB,TA));
        for e=1:numel(ib)
            yv=TB(ib(e)); x=TA(ja(e)); if x==yv, continue; end
            z=find(RB(ib(e),:) & RA(ja(e),:)); z(z==x|z==yv)=[]; hc(z)=hc(z)+1;
        end
        np=np+1;
    end
    if np>0, hc=hc/np; end
end

function m = viridis_local()
    a=[68 1 84;72 40 120;62 74 137;49 104 142;38 130 142;31 158 137; ...
       53 183 121;110 206 88;181 222 43;253 231 37]/255;
    xi=linspace(0,1,size(a,1)); xo=linspace(0,1,256);
    m=[interp1(xi,a(:,1),xo)' interp1(xi,a(:,2),xo)' interp1(xi,a(:,3),xo)'];
end

% ===================== local: ablation panel =========================
function ablation_panel(ax, data, baseline, ylo, yhi, ttl, ylab, cL, cC)
    hold(ax,'on'); off=0.16; series={-off,1,cL,'o'; off,2,cC,'s'}; hleg=[];
    for r=1:2
        d=data{series{r,2}}; c=series{r,3}; mk=series{r,4};
        for s=1:4
            [mm,lo,hi]=ci95_(d(:,s)); x0=s+series{r,1};
            plot(ax,[x0 x0],[lo hi],'-','Color',c,'LineWidth',1.4);
            plot(ax,x0+[-0.06 0.06],[lo lo],'-','Color',c,'LineWidth',1.1);
            plot(ax,x0+[-0.06 0.06],[hi hi],'-','Color',c,'LineWidth',1.1);
            hh=plot(ax,x0,mm,mk,'MarkerFaceColor',c,'MarkerEdgeColor','w','MarkerSize',6.5,'LineWidth',1.1);
            if s==1, hleg(end+1)=hh; end
        end
    end
    plot(ax,[0.5 4.5],[baseline baseline],'--','Color',[0.55 0.55 0.52],'LineWidth',1);
    set(ax,'XTick',1:4,'XTickLabel',{'drug','cell','omics','omics+net'}); xlim(ax,[0.5 4.5]); ylim(ax,[ylo yhi]);
    ylabel(ax,ylab); title(ax,ttl);
    legend(hleg,{'leave-drug-out','leave-cell-line-out'},'Location','southeast','Box','off','FontSize',8);
    fig_style(ax);
end

function [m,lo,hi]=ci95_(x)
    x=x(~isnan(x)); nn=numel(x); m=mean(x);
    if nn<2, lo=m; hi=m; return; end
    tb=[12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228];
    df=nn-1; if df>=1&&df<=10, tc=tb(df); else, tc=1.96; end
    h=tc*std(x)/sqrt(nn); lo=m-h; hi=m+h;
end

% ===================== local: plotting helpers =======================
function scatterconf(ax, x, y, cls, c0, c1)
% scatter of FFL count vs target-set size, coloured by class, with a binned
% median trend to show FFLs rise with the number of targets (the confound)
    x=x(:); y=y(:); cls=cls(:); hold(ax,'on');
    for g=0:1
        idx=find(cls==g); c=c0; if g==1, c=c1; end
        if numel(idx)>500, idx=idx(randperm(numel(idx),500)); end
        jx=(rand(numel(idx),1)-0.5)*0.8;
        try, scatter(ax, x(idx)+jx, y(idx), 9, c, 'filled','MarkerFaceAlpha',0.25,'MarkerEdgeColor','none');
        catch, plot(ax, x(idx)+jx, y(idx),'o','MarkerSize',3,'MarkerFaceColor',min(c+0.3,1),'MarkerEdgeColor','none'); end
    end
    hi=pctl(x,98); ed=linspace(min(x),hi,10); bx=[]; by=[];
    for b=1:numel(ed)-1
        m=x>=ed(b)&x<ed(b+1); if sum(m)>=5, bx(end+1)=mean(x(m)); by(end+1)=median(y(m)); end
    end
    if numel(bx)>=2, plot(ax,bx,by,'-','Color',[0.25 0.25 0.24],'LineWidth',1.8); end
    if hi>0, xlim(ax,[0 hi]); end; yh=pctl(y,97); if yh>0, ylim(ax,[0 yh]); end
end

function fig_style(ax)
    if nargin<1, ax=gca; end
    AX=[0.62 0.62 0.60]; TICKC=[0.20 0.20 0.19]; LAB=[0.05 0.05 0.05];
    set(ax,'FontName','Helvetica','FontSize',12,'LineWidth',0.75,'Box','off','TickDir','out', ...
        'TickLength',[0.012 0.012],'Layer','top','Color','w','XColor',AX,'YColor',AX);
    try, ax.XAxis.Label.Color=LAB; ax.YAxis.Label.Color=LAB; catch, end
    set(get(ax,'Title'),'FontSize',14,'FontWeight','bold','Color',LAB);
    set(get(ax,'XLabel'),'FontSize',12.5,'Color',LAB); set(get(ax,'YLabel'),'FontSize',12.5,'Color',LAB);
    try, ax.XAxis.TickLabelColor=TICKC; ax.YAxis.TickLabelColor=TICKC; catch, end
end

function [fd,xi]=kde1(z, lo, hi, ng)
    if nargin<4, ng=400; end
    xi=linspace(lo,hi,ng); z=z(:); z=z(~isnan(z));
    if exist('ksdensity','file')==2, [fd,xi]=ksdensity(z,xi); return; end
    h=1.06*std(z)*numel(z)^(-1/5); if h<=0, h=1; end
    zz=z; if numel(zz)>20000, zz=zz(randperm(numel(zz),20000)); end
    fd=zeros(1,ng); for k=1:ng, fd(k)=mean(exp(-0.5*((xi(k)-zz)/h).^2))/(h*sqrt(2*pi)); end
end

function q=pctl(x, p)
    x=sort(x(:)); n=numel(x); if n==0, q=NaN; return; end; if n==1, q=x(1); return; end
    idx=(p/100)*(n-1)+1; lo=floor(idx); hi=ceil(idx); ff=idx-lo; q=x(lo)*(1-ff)+x(hi)*ff;
end

function c=bluemap(n)
    light=[0.82 0.89 0.97]; dark=[0.09 0.28 0.55]; t=linspace(1,0,max(1,n))'; c=light+t.*(dark-light);
end

function s=addcommas(n)
    s=num2str(round(n)); neg=~isempty(s)&&s(1)=='-'; if neg, s=s(2:end); end
    out=''; c=0; for i=numel(s):-1:1, out=[s(i) out]; c=c+1; if mod(c,3)==0&&i>1, out=[',' out]; end, end
    if neg, out=['-' out]; end; s=out;
end

function boxpts(ax, vals, y, c0, c1, ttl, ylab)
    if nargin<7, ylab='value'; end
    axes(ax); hold(ax,'on'); vals=vals(:); y=y(:); m=~isnan(vals)&~isnan(y); vals=vals(m); y=y(m);
    if isempty(vals), title(ax,ttl); return; end
    cols={c0,c1}; labs={'non','syn'};
    for g=0:1
        v=vals(y==g); if isempty(v), continue; end
        x0=g+1; c=cols{g+1};
        q1=pctl(v,25); q2=pctl(v,50); q3=pctl(v,75); iqr=q3-q1; wlo=max(min(v),q1-1.5*iqr); whi=min(max(v),q3+1.5*iqr); w=0.30;
        plot(ax,[x0 x0],[q1 wlo],'-','Color',c,'LineWidth',1.1); plot(ax,[x0 x0],[q3 whi],'-','Color',c,'LineWidth',1.1);
        plot(ax,x0+[-0.10 0.10],[wlo wlo],'-','Color',c,'LineWidth',1.1); plot(ax,x0+[-0.10 0.10],[whi whi],'-','Color',c,'LineWidth',1.1);
        nn=numel(v); cap=350; if nn>cap, sel=randperm(nn,cap); vv=v(sel); else, vv=v(:); end
        jit=(rand(numel(vv),1)-0.5)*0.30; drawpts(ax, x0+jit, vv, c);
        patch(ax, x0+[-w -w w w], [q1 q3 q3 q1], c, 'FaceAlpha',0.16,'EdgeColor',c,'LineWidth',1.4);
        plot(ax, x0+[-w w],[q2 q2],'-','Color',c,'LineWidth',2.4);
        plot(ax, x0, mean(v),'d','MarkerFaceColor','w','MarkerEdgeColor',c,'MarkerSize',6.5,'LineWidth',1.2);
    end
    set(ax,'XTick',[1 2],'XTickLabel',labs); xlim(ax,[0.5 2.5]);
    lo0=min(0,pctl(vals,1)); up=pctl(vals,97); if up<=lo0, up=lo0+1; end
    rng=up-lo0; if rng<=0, rng=1; end; ylim(ax,[lo0-0.06*rng, up+0.05*rng]);
    ylabel(ax,ylab); title(ax,ttl); fig_style(ax);
end

function drawpts(ax, xx, yy, c)
    try, scatter(ax, xx, yy, 10, c, 'filled', 'MarkerFaceAlpha',0.22,'MarkerEdgeColor','none');
    catch, plot(ax, xx, yy, 'o', 'MarkerSize',3, 'MarkerFaceColor',min(c+0.30,1),'MarkerEdgeColor','none'); end
end

function [cellNames, geneNames, expr]=read_cellminer(path)
    TAB=char(9); META=6; fid=fopen(path,'r'); if fid<0, error('Cannot open %s', path); end
    cellNames={}; ncell=0;
    while true
        ln=fgetl(fid); if ~ischar(ln), fclose(fid); error('No header'); end
        f=ssplit(ln,TAB);
        if ~isempty(f)&&strncmpi(strtrim(f{1}),'Gene name',9)
            last=numel(f); while last>META&&isempty(strtrim(f{last})), last=last-1; end
            cellNames=strtrim(f(META+1:last)); ncell=numel(cellNames); break;
        end
    end
    geneNames={}; rows={};
    while true
        ln=fgetl(fid); if ~ischar(ln), break; end; if isempty(strtrim(ln)), continue; end
        f=ssplit(ln,TAB); if numel(f)<META+1, continue; end
        g=strtrim(f{1}); if isempty(g), continue; end
        v=nan(1,ncell); for c=1:ncell, idx=META+c; if idx<=numel(f), v(c)=str2double(f{idx}); end, end
        geneNames{end+1,1}=g; rows{end+1,1}=v;
    end
    fclose(fid); expr=cell2mat(rows);
end

function parts=ssplit(ln, d)
    pos=[0, find(ln==d), numel(ln)+1]; parts=cell(1,numel(pos)-1);
    for k=1:numel(pos)-1, parts{k}=ln(pos(k)+1:pos(k+1)-1); end
end

function c=normcell(s)
    c=upper(strtrim(s)); k=strfind(c,':'); if ~isempty(k), c=c(k(1)+1:end); end
    c=regexprep(c,'\(.*?\)',''); c=regexprep(c,'/ATCC',''); c=regexprep(c,'[^A-Z0-9]','');
end

function T=tg_in_signor(drug2tg, sig2idx, d)
    T=[]; if ~isKey(drug2tg,d), return; end
    g=drug2tg(d); for k=1:numel(g), if isKey(sig2idx,g{k}), T(end+1)=sig2idx(g{k}); end, end
    T=unique(T);
end

function savepdf(f, stem)
    drawnow; pdf=[stem '.pdf'];
    try, exportgraphics(f, pdf, 'ContentType','vector');
    catch
        try, set(f,'PaperPositionMode','auto','Units','inches'); p=get(f,'Position');
             set(f,'PaperUnits','inches','PaperSize',[p(3) p(4)]); print(f, pdf, '-dpdf','-r300');
        catch, print(f,[stem '.png'],'-dpng','-r300'); end
    end
    try, savefig(f,[stem '.fig']); catch, end
end

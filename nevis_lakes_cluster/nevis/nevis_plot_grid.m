function nevis_plot_grid(gg,ni)
% plot grid defined in struct gg
% if ni included, highlighs nodes labelled ni
%
% IJH 31 July 2014

if nargin<2, ni = []; end

nx = gg.nx;
ny = gg.ny;
ex = gg.ex;
ey = gg.ey;
fx = gg.fx;
fy = gg.fy;
cx = gg.cx;
cy = gg.cy;

nout = gg.nout;
eout = gg.eout;
fout = gg.fout;
cout = gg.cout;
nbdy = gg.nbdy;
ebdy = gg.ebdy;
fbdy = gg.fbdy;
cbdy = gg.cbdy;
nin = gg.nin;
ein = gg.ein;
fin = gg.fin;
cin = gg.cin;

figure;
    plot(nx(nout),ny(nout),'.','color',.9*[1 1 1],'markersize',8); hold on;
    plot(ex(eout),ey(eout),'x','color',.9*[1 1 1]);
    plot(fx(fout),fy(fout),'x','color',.9*[1 1 1]);
    plot(cx(cout),cy(cout),'+','color',.9*[1 1 1]);
    plot(nx(nbdy),ny(nbdy),'r.','markersize',8); 
    plot(nx(nin),ny(nin),'k.','markersize',8);
    plot(ex(ebdy),ey(ebdy),'rx');
    plot(ex(ein),ey(ein),'kx');
    plot(fx(fbdy),fy(fbdy),'rx');
    plot(fx(fin),fy(fin),'kx');
    plot(cx(cbdy),cy(cbdy),'r+');
    plot(cx(cin),cy(cin),'k+');
    
    plot(nx(ni),ny(ni),'go');
end
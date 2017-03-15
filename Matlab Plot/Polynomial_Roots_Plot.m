%This code generates 10000 polynomials of degree n = 100, taking into account only the polynomials that have roots smaller than 2 (which is the most of them). These roots are computed in parallel and saved in the vector z. After that, the program prints the proportion of real roots and plot the roots with density in the complex plane. This plot is also saved in the computer automatically.

z = zeros(0);
n = 101;
matlabpool('open',4);
parfor j=1:10000
   p = random('Uniform', -1, 1, [1,n]);
   R = roots(p);
   if max(abs(real(R))) < 2 & max(abs(imag(R))) < 2
      z = [ z, R.' ];
   end
end
matlabpool('close');

Re = real(z);
Im = imag(z);

c=0;
for j=1:length(Im)
   if Im(j)==0
      c=c+1;
   end
end
c/length(Im)

[values, centers] = hist3([Im(:) Re(:)],[1000 1000]);
imagesc(min(Re):.001:max(Re),min(Im):.001:max(Im), values,[0,10]);
colorbar
axis equal
axis xy
cmap = summer(max(values(:)));
cmap(1:1,:) = 0;
colormap(cmap);
set(gcf,'PaperUnits','inches','PaperSize',[1,1],'PaperPosition',[0 0 12 12]);
print('-dpng','-r100','Uniform');

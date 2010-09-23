clear all;
close all;
!cp tree.nwk tmp.nwk
!perl -i.bk -wpe's/\)([0-9]+):/,$1:100):/g' tmp.nwk
figure;
set(gcf,'Renderer','ZBuffer');
hold on;
csv=csvread('raw.csv');
csvP=exp(csvread('prior.csv'))*sum(sum(csv));
for i=1:size(csvP,1)
    for j=1:size(csvP,2)
        if csvP(i,j)<=5 && csv(i,j)<=5
            csv(i,j)=0;
        else
            csv(i,j)=log10(csv(i,j)/csvP(i,j));
        end
    end
end
csv=min(csv, 0.3);
csv=max(csv,-0.3);
csv=csv';

tree=phytreeread('tmp.nwk');
set(gcf,'Color',[1 1 1]);
subplot('Position',[0.30 0.73 0.705 0.27]);
plotX(tree,'Orientation','bottom');
subplot('Position',[0 0.00 0.3 0.705]);
[handles,order]=plotX(tree);
%csv=csv(end:-1:1,:);
im=ones(size(csv,1),size(csv,2),3);
order=order(end:-1:1);
for i=1:size(csv,1)
    for j=1:size(csv,2)
        gain=min(1,max(0, csv(str2num(order{i})+1,str2num(order{j})+1)/max(max(abs(csv)))));
        loss=min(1,max(0,-csv(str2num(order{i})+1,str2num(order{j})+1)/max(max(abs(csv)))));
        if csv(str2num(order{i})+1,str2num(order{j})+1)==0,im(size(csv,1)-i+1,j,1:3)=[0.96 0.96 0.96];continue;end
        im(size(csv,1)-i+1,j,1:3)=colorize(gain,loss);
    end
end
subplot('Position',[0.3 0 0.73 0.73]);
hold on;
axis off;
box off;
image(im);

subplot('Position',[0.05 0.8 0.2 0.1]);
for i=-100:100
    gain=min(1,max(0, i/100));
    loss=min(1,max(0,-i/100));
    im2(1,i+101,1:3)=colorize(gain,loss);
end
image(im2);
set(gca,'YColor','w');
set(gca,'XTick',[1 51 101 151 201]);
m=max(max(abs(csv)));
lab{1}=sprintf('%.2f',10^(-m));
lab{2}=sprintf('%.2f',10^(-m/2));
lab{3}='1';
lab{4}=sprintf('%.2f',10^(m/2));
lab{5}=sprintf('%.2f',10^(m));
set(gca,'XTickLabel',lab);
set(gcf,'Position',[100 100 700 700]);
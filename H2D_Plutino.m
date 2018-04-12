clear;
Filename='~/Desktop/PlutinoList.txt';
disp(['Plutino info file: ',Filename]);
fid=fopen(Filename,'r');
plinfo=textscan(fid,'%s %f %f %f %f %f %f %f %f %f');
fclose(fid);

HFilename='~/Desktop/TNOH.txt';
disp(['Plutino H file: ',HFilename]);
fid=fopen(HFilename,'r');
Hinfo=textscan(fid,'%s %f');
fclose(fid);

H2DFilename='~/Desktop/H2D.txt';
disp(['Plutino H file: ',H2DFilename]);
fid=fopen(H2DFilename,'r');
H2Dinfo=textscan(fid,'%s %f');
fclose(fid);

H2D=zeros(length(H2Dinfo{1}),1);
for i=1:length(H2Dinfo{1})
    H2D(i)=str2num(H2Dinfo{1}{i});
end
H2D=[H2D H2Dinfo{2}];
% Pl=[info{2:end}];%%% TNO Data
% H=[Hinfo{2:end}];%%% absolute magnitude
% H2D=[H2Dinfo{1:end}];%% H2D reference
Hpl=zeros(length(plinfo{1}),1);
for i=1:length(plinfo{1})
    id=find(strcmp(Hinfo{:,1},plinfo{1}{i}));
    Hpl(i)=Hinfo{2}(id);
end
Dpl=zeros(length(plinfo{1}),1);
Mpl=zeros(length(plinfo{1}),1);
Rpluto=1187;
%Mpluto=6.5607561e-9;
Mpluto=1.303e22; %kg
for i=1:length(Hpl)

    %% known plutinos
    if strcmp(plinfo{1}{i},'2004DW') % orcus
        Dpl(i)=917;
        Mpl(i)=6.41e20/Mpluto;
        disp('Orcus');
    elseif strcmp(plinfo{1}{i},'2001KX76') %ixion
        Dpl(i)=617;
        Mpl(i)=3.0e20/Mpluto;
        disp('Ixion');
    elseif strcmp(plinfo{1}{i},'2003AZ84')
        Dpl(i)=727;
        Mpl(i)=3.0e20/Mpluto;
        disp('2003AZ84');
    elseif strcmp(plinfo{1}{i},'2003VS2')
        Dpl(i)=523;
        Mpl(i)=1.5e20/Mpluto;
        disp('2003VS2');
    else
        top=find(H2D(:,1)>=Hpl(i),1,'first');
        bot=find(H2D(:,1)<=Hpl(i),1,'last');
        if top==bot
             id=top;
        else
             id=top;
        end
        Dpl(i)=H2D(id,2);
        Mpl(i)=(Dpl(i)/2/Rpluto)^3;%*Mpluto;
    end

end
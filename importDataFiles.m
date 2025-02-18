close all
clear all


%x y vx vy xi
% pathToSPHCode         = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\';
% pathToDataFileStorage = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\standardWaveTank1\dataFiles\';
% pathToImageStorage    = '\\Client\H$\ENME489ISpring2012\codes\SPH2DCPPCuda\standardWaveTank1\imageFiles\';
% pathToFffmpeg         = [pathToSPHCode 'ffmpeg01/bin/'];
% pathToHereFromFfmpeg  = ['\\Client\H$\ENME489ISpring2012\codes\DEM2DCPPCuda\standardWaveTank1\imageFiles\'];
pathToHere            = 'D:/DiabloCanyon_WCSPH_Studies/x=-3.0/';
pathToSPHCode         = 'D:/DiabloCanyon_WCSPH_Studies/x=-3.0/';
pathToDataFileStorage = [pathToSPHCode 'dataFiles/'];
pathToImageStorage    = [pathToHere 'imageFiles/'];
pathToFffmpeg         = [pathToHere 'ffmpeg01/bin/'];




storageStride  = 1;
dt             =1.0*10^(-5);


exportFrames = 1;

A = dir([pathToDataFileStorage '*.txt']);
nFiles = length(A);
%ind1 = 1;
ind1 = 1;
%get color scaling
plottingAxis = [-1.6    20   -0.0626    3.5];
%plottingAxis = [0.4067    0.5592   -0.0057    0.1118];
%plottingAxis = [ 0        1200        1300        2400];
%plottingAxis = [ 0        8800        0        2400];

plottingParams.axis = plottingAxis;

%for ind1 = [1:1:10]  
if 0
while ind1 <= nFiles

    %refresh nFiles
    A = dir([pathToDataFileStorage '*.txt']);
    nFiles = length(A);

    
   ind1 %loop iteration
   s = readInDataFile([pathToDataFileStorage A(ind1).name]);
   
   
   
   h = plot02(s,plottingAxis);
    grid on
    grid minor
    %title(['t = ' sprintf('%f',(ind1-1)*dt*storageStride) 's'],'Fontweight','normal');
    tStep = ind1*storageStride;
    %xlabel('distance (m)')
    %ylabel('height (m)')
    %niceFig()
    %plot(s(7373:end,5),'k.')
  %  plot(s(1:end,5),'k.')

   if exportFrames
   set(gcf,'visible','off','renderer','painters');
   axis equal
   xlim([-1.6    20])
   ylim([-0.0626    3.5])
   end
   
   %hold on
   %quiver(s(:,1),s(:,2),s(:,3),s(:,4),'b.')
   %set(gca,'dataaspectratio',[1 1 1]);
   %axis([-11   11   -2   15]); - original Funnel axis
   %axis([plottingAxis]); % zoom on opening
   
%   colormap copper
%colormap hot
%colormap winter
   C = colormap(jet); %define map
   C = C./max(C(:));
   colormap(C);


   caxis([0 30])
   
   hold off
     
   drawnow
   match =["dataFile0000","dataFile000","dataFile00","dataFile0"];
   str = erase(A(ind1).name,match);
   
   if exportFrames;
       cd(pathToImageStorage);
       print(gcf,'-dtiff',[str(1:end-4),'.tiff',]);
%        savefig(gcf,[str(1:end-4),'.fig',]);
%         cdata=print(gcf,'-RGBImage','-r265');
% %         image=imshow(cdata)
% %         imwrite(image,[str(1:end-4),'.tiff']);
%         F=im2frame(cdata);
%         print(F,'-dtiff',[str(1:end-4),'.tiff',]);
   end
   display([num2str(ind1) ' / ' num2str(nFiles)]);
   ind1 = ind1+1; %this does not affect "for" loop
   %pause
   %keyboard
   if 0
   cd('D:\SPH2D_Cuda\cases\Dam_Break_2D\postProcessing\')
   save(['v004timeStep'],'s')
   return
   end
   
   
end
end

if 1
    ind1=1;
    while ind1 <= nFiles
        A = dir([pathToDataFileStorage '*.txt']);
        nFiles = length(A);

    
        ind1 %loop iteration
        s = readInDataFile([pathToDataFileStorage A(ind1).name]);
   
        fprintf('saving data\n')
        I = find(s(:,4)==7);
        Freeparticles = s(I,1:2);
        save D:/DiabloCanyon_WCSPH_Studies/x=-3.0/PosFreePart4 Freeparticles
    
        DensityFreeParticles = s(I,3);
        save D:/DiabloCanyon_WCSPH_Studies/x=-3.0/DensityFreePart4 DensityFreeParticles
    
    
        I1 = find(s(:,4)==2);
        Contparticles = s(I1,1:2);
        save D:/DiabloCanyon_WCSPH_Studies/x=-3.0/PosContPart4 Contparticles
    
        DensityContainerParticles = s(I1,3);
        save D:/DiabloCanyon_WCSPH_Studies/x=-3.0/DensityContPart4 DensityContainerParticles
        
%         I2 = find(s(:,7)==5);
%         DensityMeasuredParticles = s(I2,5);
%         MeasParticles = s(I2,1:2);
%         save ..\DensityMeasPart4 DensityMeasuredParticles
        
        ind1=ind1+1;
    end
end



%path to here

% I = regexp(pathToFffmpeg,'/'); %replace / with \
% pathToFffmpeg(I) = '\'; 
% 
% I = regexp(pathToImageStorage,'/'); %replace / with \
% pathToImageStorage(I) = '\'; 
% 
%    makeMovieString = [pathToFffmpeg 'ffmpeg -i ' pathToImageStorage '\dataFile%05d.png -r 30 ' pathToImageStorage '\movieZoom1.mp4'];
%    dos(makeMovieString);
% % % % 
%    makeMovieString2 = [pathToFffmpeg 'ffmpeg -i ' pathToImageStorage '\dataFile%05d.png -r 30 -s 800x600 ' pathToImageStorage '\movieSmallZoom1.mp4'];
%    dos(makeMovieString2);
% % 
% 
% 






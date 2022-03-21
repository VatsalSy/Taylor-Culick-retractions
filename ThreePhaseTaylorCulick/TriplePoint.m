%% Tracking Three Phase Contact line (apparent)
% Last Update: November 15, 2020
% Author: Vatsal Sanjay
% vatsalsanjay@gmail.com
% Physics of Fluids
clc
clear
close
%%
ci = 3002;
plottingBool = 0;
DistCutoff = 5e-3;
dt = 5e-1;
nGFS = 1000;
fileWrite = sprintf("TriplePointTracking_%d_%3.2e.dat", ci, dt);
folder1 = sprintf('TriplePoint_%3.2e', dt); % output folder
%% Output Folder
if plottingBool == 1
  opFolder = fullfile(cd, folder1);
  if ~exist(opFolder, 'dir')
  mkdir(opFolder);
  end
end
%%
ZoomWindowR = 0.5;
rOld = 0.0;
plt.rminOld = 0.0; plt.rmaxOld = 8.0; plt.zmin = -1.0; plt.zmax = 1.0;
twrite = zeros(nGFS+1, 1); z0 = zeros(nGFS+1, 1); r0 = zeros(nGFS+1, 1);
TPcounter = 0;
%%
for ti = 1:1:nGFS
    tic
    t = (ti-1)*dt;
    filename = sprintf('intermediate/snapshot-%5.4f',t);
    if plottingBool == 1
      plt.name1 = [folder1 '/' sprintf('%4.4d.png',ti)];
    end

    if (exist(filename, 'file'))
        % Facets
        ll=evalc(sprintf('!./getFacet1 %s', filename));
        bolo=textscan(ll,'%f %f\n');
        r = bolo{2}; z = -bolo{1};
        Xf.r1 = reshape(r, [2, int32(length(r)/2)]);
        Xf.z1 = reshape(z, [2, int32(length(z)/2)]);
        ll=evalc(sprintf('!./getFacet2 %s', filename));
        bolo=textscan(ll,'%f %f\n');
        r = bolo{2}; z = -bolo{1};
        Xf.r2 = reshape(r, [2, int32(length(r)/2)]);
        Xf.z2 = reshape(z, [2, int32(length(z)/2)]);

        ll=evalc(sprintf('!./getX0Y0 %s %f %f', filename,...
            rOld, DistCutoff));
        plt.Zoom = 0;
        if ~isempty(ll)
            bolo=textscan(ll,'%f %f %f\n');
            Xf.z0 = -sum(bolo{1}./bolo{3})/sum(1./bolo{3});
            Xf.r0 = sum(bolo{2}./bolo{3})/sum(1./bolo{3});
            plt.r1 = Xf.r0 - ZoomWindowR/2; plt.r2 = Xf.r0 + ZoomWindowR/2;
            plt.z1 = Xf.z0 - ZoomWindowR/16; plt.z2 = Xf.z0 + ZoomWindowR/16;
            if (plt.r1 < 0)
                plt.r1 = 0; plt.r2 = ZoomWindowR;
            end
            plt.rmin = Xf.r0 - 4.0; plt.rmax = Xf.r0 + 4.0;
            if (plt.rmin < 0)
                plt.rmin = 0; plt.rmax = 8.0;
            end
            plt.rminOld = plt.rmin; plt.rmaxOld = plt.rmax;
            plt.Zoom = 1;
            TPcounter = TPcounter+1;
            twrite(TPcounter) = t; z0(TPcounter) = Xf.z0;
            r0(TPcounter) = Xf.r0;
            rOld = r0(TPcounter);
        else
            fprintf("Counld not get TP for %f\n", t);
            plt.r1 = -ZoomWindowR/2; plt.r2 = ZoomWindowR/2;
            plt.z1 = -ZoomWindowR/16; plt.z2 = ZoomWindowR/16;
            plt.rmin = plt.rminOld; plt.rmax = plt.rmaxOld;
        end
        fprintf('Plotting Now %d of %d\n',ti, nGFS);
        if plottingBool == 1
            plottingTime(Xf, plt);
        end
        fprintf('%d of %d\n',ti, nGFS+1);
        toc;
    else
        fprintf('%s not found. Finishing the code\n',...
            filename);
        fprintf('%d of %d\n',ti, nGFS+1);
        toc;
        break;
    end
end
twrite = twrite(1:TPcounter); r0 = r0(1:TPcounter); z0 = z0(1:TPcounter);
T = table(twrite, r0, z0);
writetable(T, fileWrite, 'Delimiter',' ')

%%
function plottingTime(Xf, plt)
ColorFacet1 = [0 204 164]/255.;
ColorFacet2 = [0.9100 0.4100 0.1700];
figure1 = figure('visible','off','WindowState','fullscreen',...
        'Color',[1 1 1]);
h1 = subplot(2,1,1);
hold on;
plot(Xf.r1, Xf.z1,'-','color',ColorFacet1,...
    'MarkerSize',30,'LineWidth',3);
plot(-Xf.r1, Xf.z1,'-','color',ColorFacet1,...
    'MarkerSize',30,'LineWidth',3);
plot(-Xf.r2, Xf.z2,'-','color',ColorFacet2,...
    'MarkerSize',30,'LineWidth',3);
plot(Xf.r2, Xf.z2,'-','color',ColorFacet2,...
    'MarkerSize',30,'LineWidth',3);
if plt.Zoom == 1
    plot([plt.r1 plt.r2], [plt.z1 plt.z1], '--','color',[0.0 0.0 0.0],...
        'MarkerSize',30,'LineWidth',3)
    plot([plt.r2 plt.r2], [plt.z1 plt.z2], '--','color',[0.0 0.0 0.0],...
        'MarkerSize',30,'LineWidth',3)
    plot([plt.r1 plt.r2], [plt.z2 plt.z2], '--','color',[0.0 0.0 0.0],...
        'MarkerSize',30,'LineWidth',3)
    plot([plt.r1 plt.r1], [plt.z1 plt.z2], '--','color',[0.0 0.0 0.0],...
        'MarkerSize',30,'LineWidth',3)
end

plot([plt.rmin plt.rmax], [plt.zmin plt.zmin], '-','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.rmax plt.rmax], [plt.zmin plt.zmax], '-','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.rmin plt.rmax], [plt.zmax plt.zmax], '-','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.rmin plt.rmin], [plt.zmin plt.zmax], '-','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)

axis equal
xlim([plt.rmin plt.rmax])
ylim([plt.zmin plt.zmax])
set(h1,'FontSize',15,'LineWidth',3);
xlabel('\boldmath{$\mathcal{R}$}','LineWidth',2,'FontWeight',...
    'bold','FontSize',20,...
    'Interpreter','latex');
ylabel('\boldmath{$\mathcal{Z}$}','LineWidth',2,'FontWeight',...
    'bold','FontSize',20,...
    'Interpreter','latex');
% axis off;

h2 = subplot(2,1,2);
hold on;
plot([plt.r1 plt.r2], [plt.z1 plt.z1], '--','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.r2 plt.r2], [plt.z1 plt.z2], '--','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.r1 plt.r2], [plt.z2 plt.z2], '--','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
plot([plt.r1 plt.r1], [plt.z1 plt.z2], '--','color',[0.0 0.0 0.0],...
    'MarkerSize',30,'LineWidth',3)
if plt.Zoom == 1
    hold on;
    plot(Xf.r1, Xf.z1,'-','color',ColorFacet1,...
        'MarkerSize',30,'LineWidth',3);
    plot(Xf.r2, Xf.z2,'-','color',ColorFacet2,...
        'MarkerSize',30,'LineWidth',3);
    plot(Xf.r0, Xf.z0, '.','color',[0.0 0.0 0.0],...
        'MarkerSize',20,'LineWidth',3)
end
axis equal
xlim([plt.r1 plt.r2])
ylim([plt.z1 plt.z2])
h2.Position = [h1.Position(1) h1.Position(2)-0.75*h1.Position(4) ...
    h2.Position(3) h2.Position(4)];
axis off;

set(figure1,'pos',[1 1 1080 960]);
export_fig(plt.name1,'-r150');
close all;
end

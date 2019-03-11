function LEiDA_EEG_plot_nodes(V)

hold on

% % PLOT CORTEX
% 
cortex.path='MNI152_T1_2mm_brain_mask.nii';
cortex.pial=mapPial(cortex.path);
cortex.color=[0.9 0.9 0.9];
cortex.transparency=0.3; % To view only opaque cortex =1;
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency); %transparency

% center origin
ori=[65 45.5 35];

load gilesCOGmni.txt gilesCOGmni
scale=5.5;
MNI_coord=scale*(gilesCOGmni/10);
clear gilesCOGmni

% Use this if you want to flip the orientation
MNI_coord(:,1)=MNI_coord(:,1)*-1;

%%  PLOT LINKS 

% Normalize the vector and assign 
% RED = community over mean
% BLUE = community below mean

thisMean=mean(V);

% if mean(V)>0
%      V=-V;
% end

n=sum(V(:)>thisMean);
if n
    [~, index] = sort(V,'descend');
    Mode2_pos=sort(index(1:n));   
    Mode1_neg=[];
    % Uncomment this to plot also links in the blue community
%     [~, index] = sort(V);  
%     Mode1_neg= sort(index(1:n));     
else
    Mode1_neg=[];
    Mode2_pos=[];
end

% Here plot red links between red nodes
for a=1:numel(Mode2_pos)
    n=Mode2_pos(a);
    for b=1:a
        p=Mode2_pos(b);
        c1=[MNI_coord(n,2)+ori(1) MNI_coord(n,1)+ori(2) MNI_coord(n,3)+ori(3)];
        c2=[MNI_coord(p,2)+ori(1) MNI_coord(p,1)+ori(2) MNI_coord(p,3)+ori(3)];
        plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r','LineWidth',1);
    end
end

% Here plot blue links between red nodes
for a=1:numel(Mode1_neg)
    n=Mode1_neg(a);
    for b=1:a
        p=Mode1_neg(b);
        c1=[MNI_coord(n,2)+ori(1) MNI_coord(n,1)+ori(2) MNI_coord(n,3)+ori(3)];
        c2=[MNI_coord(p,2)+ori(1) MNI_coord(p,1)+ori(2) MNI_coord(p,3)+ori(3)];
        plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','b','LineWidth',1);
    end
end


%% Plot shperes

% Here plot spheres all the same size but 
% scale color intensity from red to yellow or blue to cyan


[x,y,z] = sphere;
a=2; % shpere size

for n=1:length(V)
    if V(n)>0
        surf(x*a+MNI_coord(n,2)+ori(1), y*a+MNI_coord(n,1)+ori(2),z*a+MNI_coord(n,3)+ori(3),'FaceColor',[1 abs(V(n)) 0],'EdgeColor','none','FaceAlpha',1);
    elseif V(n)<0
        surf(x*a+MNI_coord(n,2)+ori(1), y*a+MNI_coord(n,1)+ori(2),z*a+MNI_coord(n,3)+ori(3),'FaceColor',[abs(V(n)) abs(V(n)) abs(V(n))],'EdgeColor','none','FaceAlpha',0.3);
    end
end

%%


% -------------------------------------------------------
% Setting image properties - light, material, angle
% -------------------------------------------------------
axis off;
axis equal
set(gcf,'Renderer', 'OpenGL') % USE UNDER LINUX FOR TRANSPARENCY 
view(3); axis off;
daspect([1 1 1]);
pbaspect([1 1 1]);
set(gca,'CameraViewAngle', 6);
set(gca, 'Projection', 'orthographic')
set(gca, 'CameraTarget', [51 68 90])
view([-90 90]) % top
material dull; lighting phong;
camlight;
rotate3d;

end

function pial=mapPial(region)

VG=spm_vol(region(1,:));
pial=zeros(VG.dim(1:3)); 
for i=1:VG.dim(3)
  pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end
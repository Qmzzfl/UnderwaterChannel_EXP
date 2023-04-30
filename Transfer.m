%%%%%%%%%% Function Definition %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name:Transfer
%Output the transfer function of the multipath channel correspond by the
%parameters
%Function parameter:
%   reflection_times:reflection times of channel
%   f_list:the 1-D matrix of argument f
%   z_source:the coordinate of source
%   z_hydrophone:the coordinate of source
%   horizon_d: horizon distance between source and hydrophone
%   Depth:the depth of sea
%   c:wave speed in media 
%   c1-media sound speed c2-seabed sound speed
%   z1-media resistance  z2-seabed resistance
%Function Return:The transfer function of the multichannel
function Transfer_Function = Transfer(reflection_times,...
                                      f_list,...
                                      z_source,z_hydrophone,...
                                      horizon_d,...
                                      Depth,c1,c2,z1,z2)
    %calculate the distance of direcet wave
    Distance_direct = sqrt((z_source-z_hydrophone)^2+horizon_d^2);
    %Define matrix of coordinate of vitual source and correspond magnitude
    z_virtualSource = zeros(2,reflection_times);
    mag_virtualSource = zeros(2,reflection_times);
    Distance = zeros(2,reflection_times);
    %Define the one-reflection parameters
    z_virtualSource(1,1) = SurfaceReflect(z_source);
    mag_virtualSource(1,1) = -1; 
    z_virtualSource(2,1) = SeabedReflect(z_source,Depth);
    mag_virtualSource(2,1) = 1; 
    % cross operation to calculate other parameters 
    for i = 2:reflection_times
        z_virtualSource(1,i) = SurfaceReflect(z_virtualSource(2,i-1));
        Distance(1,i) = sqrt((z_hydrophone-z_virtualSource(1,i)).^2+horizon_d^2);
        mag_virtualSource(1,i) = -1*mag_virtualSource(2,i-1); 

        z_virtualSource(2,i) = SeabedReflect(z_virtualSource(1,i-1),Depth);
        Distance(2,i) = sqrt((z_hydrophone-z_virtualSource(2,i)).^2+horizon_d^2);
        mag_virtualSource(2,i) = (refCo(Distance(2,i),horizon_d,z1,z2,c1,c2))^floor(i/2);
        %mag_virtualSource(2,i) = -1*mag_virtualSource(1,i-1); 
    end
    %Calculate the Distance between vitual sources and hyrdophone
    Distance = sqrt((z_hydrophone-z_virtualSource).^2+horizon_d^2);
    %Put it into 1-D for easy indexing
    Distance = reshape(Distance,1,[]);
    mag_virtualSource = reshape(mag_virtualSource,1,[]);
    
    %Accumulate to get transfer function
    Transfer_Function = 1/Distance_direct*...
            exp(-1j*2*pi*f_list*Distance_direct/c1);
    for i = 1:length(Distance)
        Transfer_Function = Transfer_Function + ...
            mag_virtualSource(i)/Distance(i)*...
            exp(-1j*2*pi*f_list*Distance(i)/c1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: SurfaceReflect
%Function parameter:the coordinate of source before reflect
%Function return:the coordinate of vitaul source after reflect by SURFACE
function z_SurfaceReflect = SurfaceReflect(z_origin)
    z_SurfaceReflect = -z_origin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: SurfaceReflect
%Function parameter:the coordinate of source before reflect
%Function return:the coordinate of vitaul source after reflect by SEABED
function z_SeabedReflect = SeabedReflect(z_origin,Depth)
    z_SeabedReflect = 2*Depth-z_origin;
end


% Function name:refCo
% Work:Calculate Reflect Coefficient in different channel
% Input:L-channel length  x-Horizontal distance  
%       z1-media resistance  z2-seabed resistance
%       c1-media sound speed c2-seabed sound speed
% Return:reflect coefficient
function Rp = refCo(L,x,z1,z2,c1,c2)
    theta_i = asin(x/L);
    theta_t = asin(c2*x/(c1*L));
    z1_ob = z1/cos(theta_i);
    z2_ob = z2/cos(theta_t);
    Rp = real((z2_ob-z1_ob)/(z2_ob+z1_ob));
end
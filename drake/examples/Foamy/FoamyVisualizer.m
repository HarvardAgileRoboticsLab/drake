classdef FoamyVisualizer < Visualizer
    % Implements the draw function for the McFoamy RC airplane
    
    properties
        
        scale = 0.3; % 30 cm length
        
        % Colors
        fuselage_color = [23 85 230]/255;
        body_color = [138 138 138]/255;
        tail_color = [126 26 26]/255;
        elev_color = [126 26 26]/255;
        wing_color = [126 26 26]/255;
        
        % Line width
        line_width = 1;
        
        % Forest
        forest = [];
        
        % Plotting limits
        xlims = [-3.5 3.5]; % meters
        ylims = [-2.5 2.5];
        zlims = [0 3];
        
        % Turf limits
        xturf = [-2 2];
        yturf = [-1.5 1.5];
        zturf = [0 .01];
        
    end
    
    methods
        function obj = FoamyVisualizer(plant)
            typecheck(plant,'FoamyPlant');
            obj = obj@Visualizer(plant.getOutputFrame);
        end
        
        function draw(obj,t,x)
            
            % Draw the plane
            persistent hFig ...
                fuselage1 fuselage2 fuselage3 fuselage4 fuselage5 ...
                body1 body2 body3 body4 body5 ...
                tail wingr wingl elevr elevl ...
                prop p len  hplane
            
            if (isempty(hFig) || ~isvalid(hFig))
                 
                hFig = sfigure(1);
                set(hFig,'DoubleBuffer', 'on');
                set(hFig,'units','normalized','outerposition',[0 0 1 1])
                
                
                x_start = 0;
                x_end = 1;
                z_start = 0.1;
                z_end = -0.1;
  
                len = abs(x_start - x_end);

                % Fuselage
                fuselage1 = obj.scale*[[x_start x_end x_end+0.1 x_end+0.1 x_end x_start]-len/2; ...
                    0 0.1   0.1   0.1 0.1 0; ...
                    z_start z_start 0.7*z_start 0.7*z_end z_end z_end];
                
                fuselage2 = obj.scale*[[x_start x_end x_end+0.1 x_end+0.1 x_end x_start]-len/2; ...
                    0 -0.1   -0.1   -0.1 -0.1 0; ...
                    z_start z_start 0.7*z_start 0.7*z_end z_end z_end];
                
                fuselage3 = obj.scale*[[x_start x_end x_end+0.1 x_end+0.1 x_end x_start]-len/2; ...
                    0 0.1   0.1   -0.1 -0.1 0; ...
                    z_start z_start 0.7*z_start  0.7*z_start z_start z_start];
                
                fuselage4 = obj.scale*[[x_start x_end x_end+0.1 x_end+0.1 x_end x_start]-len/2; ...
                    0 0.1   0.1   -0.1 -0.1 0; ...
                    z_end z_end 0.7*z_end  0.7*z_end z_end z_end];
                
                fuselage5 = obj.scale*[[x_end+0.1 x_end+0.1 x_end+0.1 x_end+0.1]-len/2; ...
                    -0.1 0.1   0.1   -0.1; ...
                    0.7*z_end 0.7*z_end -0.7*z_end  -0.7*z_end];
                
                
                % Body
                body1 = obj.scale*[[x_start+0.55*len x_start+0.9*len x_start+0.8*len x_start+0.62*len]-len/2; ...
                    0.03 0.05 0.04 0.04; ...
                    z_end z_end z_end-0.07 z_end-0.07];
                
                body2 = obj.scale*[[x_start+0.55*len x_start+0.55*len x_start+0.62*len x_start+0.62*len]-len/2; ...
                    -0.03 0.03 0.04 -0.04; ...
                    z_end z_end z_end-0.07 z_end-0.07];
                
                body3 = obj.scale*[[x_start+0.55*len x_start+0.9*len x_start+0.8*len x_start+0.62*len]-len/2; ...
                    -0.03 -0.05 -0.04 -0.04; ...
                    z_end z_end z_end-0.07 z_end-0.07];
                
                body4 = obj.scale*[[x_start+0.62*len x_start+0.8*len x_start+0.8*len x_start+0.62*len]-len/2; ...
                    0.04  0.04 -0.04 -0.04; ...
                    z_end-0.07 z_end-0.07 z_end-0.07 z_end-0.07];
                
                body5 = obj.scale*[[x_start+0.8*len x_start+0.8*len x_start+0.9*len x_start+0.9*len]-len/2; ...
                    -0.04  0.04 0.05 -0.05; ...
                    z_end-0.07 z_end-0.07 z_end z_end];
                
                
                % Tail
                tail = obj.scale*[[x_start x_start+0.2*len x_start+0.15*len x_start]-len/2; ...
                    0 0 0 0; ...
                    z_end z_end z_end-0.2 z_end-0.2];
                
                % Right Wing
                wingr = obj.scale*[[x_start+0.6*len x_start+0.85*len x_start+0.7*len x_start+0.6*len]-len/2; ...
                    0.06 0.08 0.7 0.7; ...
                    0 0 0 0];
                
                % Left Wing
                wingl = obj.scale*[[x_start+0.6*len x_start+0.85*len x_start+0.7*len x_start+0.6*len]-len/2; ...
                    -0.06 -0.08 -0.7 -0.7; ...
                    0 0 0 0];
                
                % Elevator
                elevr = obj.scale*[[x_start x_start x_start+0.1*len x_start+0.15*len]-len/2; ...
                    0 0.3 0.3 0; ...
                    0 0 0 0];
                
                elevl = obj.scale*[[x_start x_start x_start+0.1*len x_start+0.15*len]-len/2; ...
                    0 -0.3 -0.3 0; ...
                    0 0 0 0];
                
                % Propeller
                a = linspace(0,2*pi,50);
                L = 0.3;
                r = 0.3;
                Rprop = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
                prop = Rprop*[L*cos(a);.02*sin(2*(a))];
                
                p = obj.scale*[repmat(x_end+0.11,1,length(prop(1,:)))-len/2;r*prop(1,:);prop(2,:)];
                
            end
          
            
            sfigure(hFig); % cla;
            hold on, grid on, view(-210,10); % set viewing angle
            if exist('hplane','var')
                try
                delete(hplane);
                catch
                end
            end
            
            % convert quaternion to rotation matrix
            xyz = x(1:3);
            q = x(4:7);
            RotMat = eye(3) + 2*hat(q(2:4))*(hat(q(2:4)) + eye(3)*q(1));
            
            function xhat = hat(x)
            % Computes the hat map, mapping a 3-vector into a 3x3 skew-symmetric matrix
            % equivalent to the cross-product operation

            xhat = [  0   -x(3)  x(2)
                     x(3)   0   -x(1)
                    -x(2)  x(1)  0];
            end
            
            % Rotation matrix applied to 3x6 persistent object vector, then
            % translated by xyz vector (repeated to match 3x6 vector size)
            fuselage1_plot = RotMat*(fuselage1) + repmat(xyz,1,length(fuselage1));
            fuselage2_plot = RotMat*(fuselage2) + repmat(xyz,1,length(fuselage2));
            fuselage3_plot = RotMat*(fuselage3) + repmat(xyz,1,length(fuselage3));
            fuselage4_plot = RotMat*(fuselage4) + repmat(xyz,1,length(fuselage4));
            fuselage5_plot = RotMat*(fuselage5) + repmat(xyz,1,length(fuselage5));
            
            body1_plot = RotMat*(body1) + repmat(xyz,1,length(body1));
            body2_plot = RotMat*(body2) + repmat(xyz,1,length(body2));
            body3_plot = RotMat*(body3) + repmat(xyz,1,length(body3));
            body4_plot = RotMat*(body4) + repmat(xyz,1,length(body4));
            body5_plot = RotMat*(body5) + repmat(xyz,1,length(body5));
            
            tail_plot = RotMat*(tail) + repmat(xyz,1,length(tail));
            wingr_plot = RotMat*(wingr) + repmat(xyz,1,length(wingr));
            wingl_plot = RotMat*(wingl) + repmat(xyz,1,length(wingl));
            elevr_plot = RotMat*(elevr) + repmat(xyz,1,length(elevr));
            elevl_plot = RotMat*(elevl) + repmat(xyz,1,length(elevl));
            p_plot = RotMat*(p) + repmat(xyz,1,length(p));
            
            hplane = [];
            hplane = [hplane, patch(fuselage1_plot(1,:),fuselage1_plot(2,:),fuselage1_plot(3,:),obj.fuselage_color,'LineWidth',obj.line_width)];
            hplane = [hplane, patch(fuselage2_plot(1,:),fuselage2_plot(2,:),fuselage2_plot(3,:),obj.fuselage_color,'LineWidth',obj.line_width)];
            hplane = [hplane, patch(fuselage3_plot(1,:),fuselage3_plot(2,:),fuselage3_plot(3,:),obj.fuselage_color,'LineWidth',obj.line_width)];
            hplane = [hplane, patch(fuselage4_plot(1,:),fuselage4_plot(2,:),fuselage4_plot(3,:),obj.fuselage_color,'LineWidth',obj.line_width)];
            hplane = [hplane, patch(fuselage5_plot(1,:),fuselage5_plot(2,:),fuselage5_plot(3,:),obj.fuselage_color,'LineWidth',obj.line_width)];
            
            hplane = [hplane, patch(body1_plot(1,:),body1_plot(2,:),body1_plot(3,:),obj.body_color,'LineWidth',obj.line_width,'FaceAlpha',0.7)];
            hplane = [hplane, patch(body2_plot(1,:),body2_plot(2,:),body2_plot(3,:),obj.body_color,'LineWidth',obj.line_width,'FaceAlpha',0.7)];
            hplane = [hplane, patch(body3_plot(1,:),body3_plot(2,:),body3_plot(3,:),obj.body_color,'LineWidth',obj.line_width,'FaceAlpha',0.7)];
            hplane = [hplane, patch(body4_plot(1,:),body4_plot(2,:),body4_plot(3,:),obj.body_color,'LineWidth',obj.line_width,'FaceAlpha',0.7)];
            hplane = [hplane, patch(body5_plot(1,:),body5_plot(2,:),body5_plot(3,:),obj.body_color,'LineWidth',obj.line_width,'FaceAlpha',0.7)];
            
            hplane = [hplane, patch(tail_plot(1,:),tail_plot(2,:),tail_plot(3,:),obj.tail_color,'LineWidth',obj.line_width)];
            
            hplane = [hplane, patch(wingr_plot(1,:),wingr_plot(2,:),wingr_plot(3,:),obj.wing_color,'LineWidth',obj.line_width)];
            hplane = [hplane, patch(wingl_plot(1,:),wingl_plot(2,:),wingl_plot(3,:),obj.wing_color,'LineWidth',obj.line_width)];
            
            hplane = [hplane, patch(elevr_plot(1,:),elevr_plot(2,:),elevr_plot(3,:),obj.elev_color,'LineWidth',obj.line_width)];
            
            hplane = [hplane,  patch(elevl_plot(1,:),elevl_plot(2,:),elevl_plot(3,:),obj.elev_color,'LineWidth',obj.line_width)];
            
            hplane = [hplane, patch(p_plot(1,:),p_plot(2,:),p_plot(3,:),'k','LineWidth',obj.line_width)];
            
            
            % plot3(x(1),x(2),x(3),'b.') % plot trajectory if desired
            % GROUND
            patch([obj.xlims(1) obj.xlims(1) obj.xlims(2) obj.xlims(2)], ...
                [obj.ylims(1) obj.ylims(2) obj.ylims(2) obj.ylims(1)], [0 0 0 0], [0.2 0.2 0.2])
            % TURF
            patch([obj.xturf(1) obj.xturf(1) obj.xturf(2) obj.xturf(2)], ...
                [obj.yturf(1) obj.yturf(2) obj.yturf(2) obj.yturf(1)], [0.01 0.01 0.01 0.01], [0 0.65 0.25])
            % WHITE LINES
            plot3([obj.xturf(1) obj.xturf(2)], [0 0], [0.01 0.01], 'w-', 'linewidth', 4), 
            plot3([0 0], [obj.yturf(1) obj.yturf(2)], [0.01 0.01], 'w-', 'linewidth', 4)
            axis image
            axis([obj.xlims(1) obj.xlims(2) obj.ylims(1) obj.ylims(2) ...
                obj.zlims(1) obj.zlims(2)]);
            title(['t = ', num2str(t(1),'%.2f') ' sec']);
            xlabel('X [m]', 'fontsize', 14)
            ylabel('Y [m]', 'fontsize', 14)
            zlabel('Z [m]', 'fontsize', 14)
            hold on
            drawnow;
        end
    end
    
end

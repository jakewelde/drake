classdef PlanarRocketVisualizer < Visualizer
% Implements the draw function for the Planar Quadrotor model

  properties
    rocket = [];
    platform = {};
  end

  methods
    function obj = PlanarRocketVisualizer(plant)
      typecheck(plant,'PlanarRocketPlant');
      obj = obj@Visualizer(plant.getOutputFrame);
      obj.rocket = plant;
    end
    
    function obj = addPlatform(obj, x, y, width)
      obj.platform = [obj.platform struct('x', x, 'y', y, 'width', width)];
    end
    
    function draw(obj,t,x)
      % Draw the rocket!!!  
      persistent hFig base motor;

      if (isempty(hFig))
        hFig = sfigure(25);
        set(hFig,'DoubleBuffer', 'on');
        
        base = [obj.rocket.wdraw*[1 -1 -1 1]/2; obj.rocket.ldraw*[1 1 -1 -1]/2];
        motor = [obj.rocket.wdraw*1.2*[1 1 -1 -1]; obj.rocket.wdraw*[1 0 0 1]];
        a = linspace(0,2*pi,50);
      end
        
      sfigure(hFig); cla; hold on; view(0,90);
      
      r = [cos(x(3)), sin(x(3)); -sin(x(3)), cos(x(3))];
      
      p = r*base;
      patch(x(1)+p(1,:), x(2)+p(2,:),1+0*p(1,:),'b','FaceColor',[.6 .6 .6])
      
      p = r*(repmat(obj.rocket.T, [1, 4]) + [motor(1,:);motor(2,:)]);
      patch(x(1)+p(1,:),x(2)+p(2,:),0*p(1,:),'b','FaceColor',[0 0 0]);
      
      if (~isempty(obj.platform))
        for i=1:numel(obj.platform)
          pl = obj.platform{i};
          patch([pl.x-pl.width/2; pl.x+pl.width/2; pl.x+pl.width/2; pl.x-pl.width/2], ...
                [pl.y-50; pl.y-50; pl.y; pl.y], ...
                'b', 'FaceColor', [0 0 1.0]);
        end
      end
      
      title(['t = ', num2str(t(1),'%.2f') ' sec']);
      set(gca,'XTick',[],'YTick',[])
      
      axis image; axis([-40.0 40.0 -1.0 60.0]);
      drawnow;
    end    
  end
  
end

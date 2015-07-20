perturbs = {};
agreements = {};
success_miqps = {};
success_reacts = {};
step_orders_miqp = {};
step_orders_react = {};

while (1)
   perturb = 2*rand([2, 1])-1;
   perturb = perturb*0.9;
   perturbs = [perturbs; perturb];
   example_options = struct('T', 2.5);
   try
       [success_miqp, step_ordering_miqp] = runAtlasStepRecovery(perturb);
   catch
       success_miqp = 0;
       step_ordering_miqp = [];
   end
   try
       [success_react, step_ordering_react] = drakeReactiveRecovery([perturb; 0.0]);
   catch
       success_react = 0;
       step_ordering_react = [];
   end
   success_miqps = [success_miqps; success_miqp];
   step_orders_miqp = [step_orders_miqp; {step_ordering_miqp}];
   success_reacts = [success_reacts; success_react];
   step_orders_react = [step_orders_react; {step_ordering_react}];
   agreement = 0;
   if (success_miqp == success_react)
       if (success_miqp == 0)
           agreement = 1; % agreed on failure
       elseif (numel(step_ordering_miqp) == numel(step_ordering_react))
           agreement = 1;
           for i=1:numel(step_ordering_miqp)
               if (size(step_ordering_miqp{i}) ~= size(step_ordering_react{i}))
                   agreement = 0;
                   break;
               end
               if (~all(step_ordering_miqp{i} == step_ordering_react{i}))
                   agreement = 0;
                   break;
               end
           end
       end
   end
   agreements = [agreements; agreement];
   
   save(['recovery_testing_en_masse' datestr(now) '.mat'], 'agreements', 'perturbs', ...
       'step_orders_miqp', 'step_orders_react', 'success_miqps', 'success_reacts');
end


figure;
perturbs_plot = horzcat(perturbs{:}).';
success_miqps_plot = horzcat(success_miqps{:}).' == 1;
success_miqps_pts = perturbs_plot(success_miqps_plot, :);
success_reacts_plot = horzcat(success_reacts{:}).' == 1;
success_reacts_pts = perturbs_plot(success_reacts_plot, :);
success_both_plot = success_miqps_plot & success_reacts_plot;
success_both_pts = perturbs_plot(success_both_plot, :);

color_neither = [0 0 0];
colors = repmat(color_neither, [numel(perturbs), 1]);
color_both = [0 1 0];
color_miqp_only = [0 0 1];
color_reacts_only = [1 0 0];
for i=1:numel(perturbs)
    if (success_miqps_plot(i) && success_reacts_plot(i))
        colors(i, :) = color_both;
    elseif (success_miqps_plot(i))
        colors(i, :) = color_miqp_only;
    elseif (success_reacts_plot(i))
        colors(i, :) = color_reacts_only;
    end
end % black = none of them worked

scatter(perturbs_plot(:, 1), perturbs_plot(:, 2), 50, colors);
hold on;
        
% all chull for pts
chull_perturbs_miqps = convhull(perturbs_plot);
line(perturbs_plot(chull_perturbs_miqps, 1), perturbs_plot(chull_perturbs_miqps, 2), ... 
    'color',color_neither);

% successful chull for miqps
chull_successful_miqps = convhull(success_miqps_pts);
line(success_miqps_pts(chull_successful_miqps, 1), success_miqps_pts(chull_successful_miqps, 2), ...
    'color',color_miqp_only);

% successful chull for reacts
chull_successful_reacts = convhull(success_reacts_pts);
line(success_reacts_pts(chull_successful_reacts, 1), success_reacts_pts(chull_successful_reacts, 2), ...
    'color',color_reacts_only);


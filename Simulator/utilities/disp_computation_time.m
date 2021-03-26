function [hours, minutes, seconds] = disp_computation_time( computation_time )

seconds = floor(mod(computation_time, 60));
minutes = floor(mod(computation_time/60, 60));
hours   = floor(mod(computation_time/60/60, 60));
disp(['Computation time: ', ...
      num2str(hours),   'h ',       ...
      num2str(minutes), 'm ',       ...
      num2str(seconds), 's.' ]);
end
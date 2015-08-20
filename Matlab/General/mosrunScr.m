num_processes = 12; 

rand('state',sum(100*clock))
process_id = round(1e4*rand(1));
for i=1:num_processes,
  file_i = ['mosrun_', num2str(i)];
  unix(['echo "step = ', num2str(num_processes), '; start = ', num2str(i), ';" >! ./tmp/', file_i]);
  unix(['cat ./tmp/', file_i, ' check_sum_marg.m >! ./tmp/', file_i, '.m']);
  
  unix(['mosrun -q -b -E"`pawd`/tmp" -J', num2str(process_id), ...
        ' /cs/phd/cheny/scripts/matlab -nojvm -nodesktop -nodisplay -r ', ...
        file_i, ' >&! ./tmp/log_', file_i, ' &']);

  unix('sleep 20');
end

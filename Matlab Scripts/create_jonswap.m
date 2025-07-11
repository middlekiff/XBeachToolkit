
% This script prompts for JONSWAP wave parameters and writes them 
% to a text file in the required format for XBeach.

% Construct the full file path in mod_dir
% Don't change the name of the file without also changing the parameters
% file to match, otherwise the simulation will fail.

gammajsp = 3.3; % Peak enhancement factor gammajsp
s = 20; % Directional spreading coefficient
fnyq = 1; % Highest frequency fnyq


filename = fullfile(outdir, 'jonswap.txt');

% Open file for writing 
fid = fopen(filename, 'w'); 
if fid == -1 
    error('Error: Cannot open file %s for writing.', filename); 
end

% Write parameters to file in the desired format 
fprintf(fid, 'Hm0 = %g\n', Hm0); 
fprintf(fid, 'Tp = %g\n', Tp); 
fprintf(fid, 'mainang = %g\n', mainang); 
fprintf(fid, 'gammajsp = %g\n', gammajsp); 
fprintf(fid, 's = %g\n', s); 
fprintf(fid, 'fnyq = %g\n', fnyq);

% Close the file 
fclose(fid);

fprintf('JONSWAP file "%s" created successfully.\n', filename);
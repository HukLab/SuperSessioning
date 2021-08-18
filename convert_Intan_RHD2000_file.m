function recParameter=convert_Intan_RHD2000_file(RawFile,OutFile,recParameter)

fid = fopen(RawFile, 'r');

s = dir(RawFile);
filesize = s.bytes;

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');

fprintf(1, '\n');
fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
    data_file_main_version_number, data_file_secondary_version_number);
fprintf(1, '\n');

if (data_file_main_version_number == 1)
    num_samples_per_data_block = 60;
else
    num_samples_per_data_block = 128;
end

Nchx=sum(recParameter.ChMask,2);    

% Read information of sampling rate and amplifier frequency settings.
sample_rate = fread(fid, 1, 'single');
recParameter.sRate=sample_rate;
dsp_enabled = fread(fid, 1, 'int16');
actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
actual_lower_bandwidth = fread(fid, 1, 'single');
actual_upper_bandwidth = fread(fid, 1, 'single');

desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
desired_lower_bandwidth = fread(fid, 1, 'single');
desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
notch_filter_mode = fread(fid, 1, 'int16');
notch_filter_frequency = 0;
if (notch_filter_mode == 1)
    notch_filter_frequency = 50;
elseif (notch_filter_mode == 2)
    notch_filter_frequency = 60;
end

desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');

% Place notes in data strucure
notes = struct( ...
    'note1', fread_QString(fid), ...
    'note2', fread_QString(fid), ...
    'note3', fread_QString(fid) );
    
% If data file is from GUI v1.1 or later, see if temperature sensor data
% was saved.
num_temp_sensor_channels = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
    || (data_file_main_version_number > 1))
    num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load board mode.
board_mode = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
    || (data_file_main_version_number > 1))
    board_mode = fread(fid, 1, 'int16');
end

% If data file is from v2.0 or later (Intan Recording Controller),
% load name of digital reference channel.
if (data_file_main_version_number > 1)
    reference_channel = fread_QString(fid);
end

% Place frequency-related information in data structure.
frequency_parameters = struct( ...
    'amplifier_sample_rate', sample_rate, ...
    'aux_input_sample_rate', sample_rate / 4, ...
    'supply_voltage_sample_rate', sample_rate / num_samples_per_data_block, ...
    'board_adc_sample_rate', sample_rate, ...
    'board_dig_in_sample_rate', sample_rate, ...
    'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
    'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
    'dsp_enabled', dsp_enabled, ...
    'desired_lower_bandwidth', desired_lower_bandwidth, ...
    'actual_lower_bandwidth', actual_lower_bandwidth, ...
    'desired_upper_bandwidth', desired_upper_bandwidth, ...
    'actual_upper_bandwidth', actual_upper_bandwidth, ...
    'notch_filter_frequency', notch_filter_frequency, ...
    'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
    'actual_impedance_test_frequency', actual_impedance_test_frequency );

% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
    'voltage_trigger_mode', {}, ...
    'voltage_threshold', {}, ...
    'digital_trigger_channel', {}, ...
    'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
    'native_channel_name', {}, ...
    'custom_channel_name', {}, ...
    'native_order', {}, ...
    'custom_order', {}, ...
    'board_stream', {}, ...
    'chip_channel', {}, ...
    'port_name', {}, ...
    'port_prefix', {}, ...
    'port_number', {}, ...
    'electrode_impedance_magnitude', {}, ...
    'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);



% Create structure arrays for each type of data channel.
amplifier_channels = struct(channel_struct);
aux_input_channels = struct(channel_struct);
supply_voltage_channels = struct(channel_struct);
board_adc_channels = struct(channel_struct);
board_dig_in_channels = struct(channel_struct);
board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;

% Read signal summary from data file header.

number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:number_of_signal_groups
    signal_group_name = fread_QString(fid);
    signal_group_prefix = fread_QString(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
        new_channel(1).port_name = signal_group_name;
        new_channel(1).port_prefix = signal_group_prefix;
        new_channel(1).port_number = signal_group;
        for signal_channel = 1:signal_group_num_channels
            new_channel(1).native_channel_name = fread_QString(fid);
            new_channel(1).custom_channel_name = fread_QString(fid);
            new_channel(1).native_order = fread(fid, 1, 'int16');
            new_channel(1).custom_order = fread(fid, 1, 'int16');
            signal_type = fread(fid, 1, 'int16');
            channel_enabled = fread(fid, 1, 'int16');
            new_channel(1).chip_channel = fread(fid, 1, 'int16');
            new_channel(1).board_stream = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
            
            if (channel_enabled)
                switch (signal_type)
                    case 0
                        amplifier_channels(amplifier_index) = new_channel;
                        spike_triggers(amplifier_index) = new_trigger_channel;
                        amplifier_index = amplifier_index + 1;
                    case 1
                        aux_input_channels(aux_input_index) = new_channel;
                        aux_input_index = aux_input_index + 1;
                    case 2
                        supply_voltage_channels(supply_voltage_index) = new_channel;
                        supply_voltage_index = supply_voltage_index + 1;
                    case 3
                        board_adc_channels(board_adc_index) = new_channel;
                        board_adc_index = board_adc_index + 1;
                    case 4
                        board_dig_in_channels(board_dig_in_index) = new_channel;
                        board_dig_in_index = board_dig_in_index + 1;
                    case 5
                        board_dig_out_channels(board_dig_out_index) = new_channel;
                        board_dig_out_index = board_dig_out_index + 1;
                    otherwise
                        error('Unknown channel type');
                end
            end
            
        end
    end
end

% Summarize contents of data file.
num_amplifier_channels = amplifier_index - 1;
num_aux_input_channels = aux_input_index - 1;
num_supply_voltage_channels = supply_voltage_index - 1;
num_board_adc_channels = board_adc_index - 1;
num_board_dig_in_channels = board_dig_in_index - 1;
num_board_dig_out_channels = board_dig_out_index - 1;

fprintf(1, 'Found %d amplifier channel.\n', ...
    num_amplifier_channels);
fprintf(1, 'Found %d auxiliary input channel.\n', ...
    num_aux_input_channels);
fprintf(1, 'Found %d supply voltage channel.\n', ...
    num_supply_voltage_channels);
fprintf(1, 'Found %d board ADC channel.\n', ...
    num_board_adc_channels);
fprintf(1, 'Found %d board digital input channel.\n', ...
    num_board_dig_in_channels);
fprintf(1, 'Found %d board digital output channel.\n', ...
    num_board_dig_out_channels);
fprintf(1, 'Found %d temperature sensor channel.\n', ...
    num_temp_sensor_channels);
fprintf(1, '\n');

% Determine how many samples the data file contains.

% Each data block contains num_samples_per_data_block amplifier samples.
bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
% Auxiliary inputs are sampled 4x slower than amplifiers
bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
% Supply voltage is sampled once per data block
bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num_board_dig_in_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num_board_dig_out_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end
% Temp sensor is sampled once per data block
if (num_temp_sensor_channels > 0)
   bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
end

% How many data blocks remain in this file?
data_present = 0;
bytes_remaining = filesize - ftell(fid);
if (bytes_remaining > 0)
    data_present = 1;
end

num_data_blocks = bytes_remaining / bytes_per_block;

num_amplifier_samples = num_samples_per_data_block * num_data_blocks;

record_time = num_amplifier_samples / sample_rate;

if (data_present)
    fprintf(1, 'File contains %0.3f seconds of data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
        record_time, sample_rate / 1000);
    fprintf(1, '\n');
else
    fprintf(1, 'Header file contains no data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
        sample_rate / 1000);
    fprintf(1, '\n');
end

%recParameter.readfromKWIK
recParameter.bitVolt=0.195;
recParameter.HdfRawDataPath='/data';
if recParameter.tEnd==-1
    recParameter.nEnd=double(num_amplifier_samples);
else
    recParameter.nEnd=double(floor(recParameter.tEnd*60*sample_rate));
end
if recParameter.tStart==0
    recParameter.nStart=1;
else
    recParameter.nStart=double(floor(recParameter.tStart*60*sample_rate)+1);
end

if (data_present)
    
    % Pre-allocate memory for data.
    fprintf(1, 'Allocating memory for data...\n');
    %output .hdf
    h5create(OutFile,recParameter.HdfRawDataPath,[Nchx recParameter.nEnd-recParameter.nStart+1],'ChunkSize',[1 num_samples_per_data_block],'Datatype','int16','FillValue',int16(0))
    %not sure whether needed
    %h5create(OutFile,'/time',recParameter.nEnd-recParameter.nStart+1,'ChunkSize',num_samples_per_data_block,'Datatype','uint32','FillValue',int32(0))
    t_amplifier = zeros(1, recParameter.nEnd);

    % Read sampled data from file.
    fprintf(1, 'Reading data from file...\n');
    
    amplifier_index = 1;
    channel_cutout=recParameter.ChMap(recParameter.ChMask);
    nSkip=0;
    if (num_aux_input_channels > 0)
        nSkip=nSkip+((num_samples_per_data_block / 4)*num_aux_input_channels)*2;
    end
    if (num_supply_voltage_channels > 0)
        nSkip=nSkip+num_supply_voltage_channels*2;
    end
    if (num_temp_sensor_channels > 0)
        nSkip=nSkip+num_temp_sensor_channels*2;
    end
    if (num_board_adc_channels > 0)
        nSkip=nSkip+(num_samples_per_data_block*num_board_adc_channels)*2;
    end
    if (num_board_dig_in_channels > 0)
        nSkip=nSkip+num_samples_per_data_block*2;
    end
    if (num_board_dig_out_channels > 0)
        nSkip=nSkip+num_samples_per_data_block*2;
    end
    for i=1:num_data_blocks
        if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
        || (data_file_main_version_number > 1))
            t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
        else
            t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint32');
        end
        %h5write(OutFile,'/time',uint32(t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)))),...
        %                amplifier_index,num_samples_per_data_block);
        if (num_amplifier_channels > 0)
            amplifier_data= fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
            h5write(OutFile,recParameter.HdfRawDataPath,int16(amplifier_data(channel_cutout,:) - 32768),...
                        [1,amplifier_index],[Nchx,num_samples_per_data_block]);
        end
        fseek(fid,nSkip,0);
        amplifier_index = amplifier_index + num_samples_per_data_block;
    end

    % Make sure we have read exactly the right amount of data.
    bytes_remaining = filesize - ftell(fid);
    if (bytes_remaining ~= 0)
        %error('Error: End of file not reached.');
    end

end

% Close data file.
fclose(fid);


% Check for gaps in timestamps.
num_gaps = sum(diff(t_amplifier) ~= 1);
if (num_gaps == 0)
    fprintf(1, 'No missing timestamps in data.\n');
else
    fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', ...
        num_gaps);
end
return


function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return



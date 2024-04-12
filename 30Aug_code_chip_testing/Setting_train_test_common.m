tim_chk=clock;disp(['Starting config at ' num2str(round(tim_chk(5:6)))]);

saveOption = 1;
plotOption = 1;

srcPosition = savePosition;
modelName = strcat('mnist_model_DELM65n_relu_dvdd_1.2_T_50_expand_8.mat');
modelFile = strcat(srcPosition, modelName);
activation = 'relu';
ser_or_par='par'; %for beta only par working so never used ser

mclk_div = 2;  % mclk cycle period:   FPGA output freq/(2*(mclk_div+1))
sclk_div = 2;  % sck cycle period:   FPGA output freq/(2*(sck_div+1))
SEL = 0;   % unsigned 2b, 0: external; 1: DLL; 2: CCO; 3: back-up; only 0 is working

pipesize_in=128;
pipesize_out_h=256;
pipesize_out_mac=64;
pipesize_beta=1025*16;

DELM65n_global_setting;

mdac = 127;  % 8b unsigned
if (exist('shift_mdac_en')==0)
	shift_mdac_en = 0;  % 1b unsigned to decide if shift_mdac on or off; since DAC settling not dominant better switch OFF
	shift_mdac = 0;  % 7b unsigned to decide shift_mdac current
end
if (exist('shift_cancel_mdac_en')==0)
	shift_cancel_mdac_en = 0;% 1b unsigned to decide if shift_cancel_mdac on or off; since DAC settling not dominant better switch OFF; can still use for creating negative bias at hidden neuron
	shift_cancel_mdac = 0;  % 7b unsigned to decide shift_cancel_mdac current
end
if (exist('out_bankB')==0)
	out_bankB = 0;	% 1b unsigned; 0 for enable bank with sigma=0.75 as W/L=290/120; both A & B can be 1
end
if (exist('out_bankA')==0)
	out_bankA = 1; % 1b unsigned; 0 for enable bank with sigma=0.6 as W/L=360/150; both A & B can be 1
end
if (exist('out_reduce')==0)
	out_reduce = 0; 	% 1b unsigned; 0(1) for high(low) power high(low) speed;
end
mu_shift = 7; % 4b unsigned; to tell divider for calculating mu from sum of values
if (exist('mu_select')==0)
	mu_select = 0; 	% 1b unsigned
end
if (exist('mu_value')==0)
	mu_value = 0;   % 8b signed
end
mccodac_en = 0;  % 1b unsigned to decide if mccodac on or off; since hold through this not proven better switch it OFF
mccodac = 0;  % 7b unsigned to decide mccodac current

%ser_or_par='ser';	% works but don't use
DELM65n_ic_input_config; %%% Configuration phase depending on above defined parameters for digi_input_layer

if (exist('out_cnt_ana')==0)
	out_cnt_ana = 100; % 15b unsigned %out_cnt_ana=1 for all positive h, otherwise h are all negative
end
if (exist('out_cnt_stl')==0)
	out_cnt_stl = 20; % 16b unsigned
end
if (exist('out_cnt_hln')==0)
	out_cnt_hln = 20; % 16b unsigned
end
reg_g_threshold = 0; % 12b unsigned

switch activation  % load from modelFile
	case 'relu'
		reg_g_threshold = 2047;
		reg_func_g_mode = 0; % 2b unsigned, 0: linear rectifier with saturation
	case 'relu_satu'
		reg_func_g_mode = 0; % 2b unsigned, 0: linear rectifier with saturation
	case 'abs_satu'
		reg_func_g_mode = 1; % 2b unsigned, 1: absolute with saturation
	case 'abs_tristate'
		reg_func_g_mode = 2; % 2b unsigned, 2 or 3: tristate
	case 'tristate'
		reg_func_g_mode = 3; % 2b unsigned, 2 or 3: tristate
end

switch h_or_mac
	case 'h'
		reg_h_mac_puf = 1; % 2b unsigned, 0: mac; 1: h;
	case 'mac'
		reg_h_mac_puf = 0;
end

bad = [47 52]; %can use if want to switch of particular physical hidden neuron
sel_hidden = [];
for i=1:128
	if(isempty(find(bad==i, 1)))
		sel_hidden = [sel_hidden, i];
	end
end
bad=[];sel_hidden = 1:128;%override above and just select all physical hidden neuron
out_ColSel_m = 0;

out_satu = 3; % unsigned 2b, 3: 1023; 2: 511; 1: 255; 0: 127.
out_Cap = 0; % unsigned 2b
out_E_cap = 0; % unsigned 1b
out_E_RX2 = 0; % unsigned 1b
if (exist('out_OP')==0)
	out_OP = 2; % unsigned 2b
end
ctrl_hold = 0; % unsigned 1b
out_Cap_m = 0; % unsigned 2b
out_satu_m = 3; % unsigned 2b
out_cfg_DLL = 6;    % for DLL, 15: X1; 6: X4; 10: X6; 12: X8
% for CCO, 127: ~125MHz; 63: ~75MHz; 31: ~40MHz;
if(SEL)
	reg_spi_clk_speed = 3; % usigned 4b, odd number only: clk_syn = mclk/(reg_spi_clk_speed+1)
else
	reg_spi_clk_speed = 0; % usigned 4b, odd number only: clk_syn = mclk/(reg_spi_clk_speed+1)
end
REG_class = N_class; % unsigned 4 bit, number of output neurons
REG_dist = 1023;   % unsigned 10 bit, default: 1023,

reg_puf_thres = 511;  % 12b
chal=uint8(round(rand(1,1)*127)); % unsigned 7b
seed=uint8(round(rand(1,1)*127)); % unsigned 7b

%ser_or_par='ser';	% works but don't use
DELM65n_ic_output_config; %%% Configuration phase depending on above defined parameters for digi_output_layer
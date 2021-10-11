function savePth = ogDownloadData(savePth, DataToDownload, ConditionName, whichDataTypes)
%
% Download data from the paper:
%   Eline R. Kupers, Marisa Carrasco, Jonathan Winawer. (2019) 
% Modeling visual performance differences 'around' the visual field: A computational
% observer approach (PLOS Computational Biology. VOLUME. ISSUE. DOI.)
%
% Data can either be cone absorptions or classification accuracy. This
% function will automatically download the latest version for each
% dataset from the Open Science Framework website.
% 
% savePth = ogDownloadData(savePth, DataToDownload, ConditionName)

% ----------------------------------------------------------------
% INPUTS:
% -----------------
%   savePth: Path to store data.
%                   [default = fullfile(ogRootPath,'data')];
%
%   DataToDownload: String that specifies whether to download absorptions
%                   or classification accuracy.
%                   Choose: 'classification' or 'absorptions'
%
%   ConditionName:  String that speficies what experimental condition to
%                   download data from.
%                   Choose: 'conedensity','conetypes', 'conetypesmixed', 
%                   'defaultnophaseshift', 'defaultnophaseshift-800trials',
%                   'defocus', 'eyemov', 'eyemovnophaseshift', 'idealobserver' 
%
%   whichDataTypes: Cell array of one or more of types of datasets.
%                   choose: {'average' 'run1' 'run2' 'run3', 'run4',
%                   'run5'}
%                   
%                   [default='average']
%
% OUTPUTS:
% -----------------
%   savePth: path where data was written
%
% EXAMPLES:
% -----------------
% Example 1: Download average classification data from defocus experiment
%   savePth = ogDownloadData('classification', 'defocus', {'average'});
% Example 2: Download raw data from all subjects
%   savePth = ogDownloadData('absorptions', 'conetypesmixed', {'run1'});

% Argument check
if notDefined('savePth'),        savePth = fullfile(ogRootPath, 'data'); end
if notDefined('whichDataTypes'),  whichDataTypes = {'average'}; end

% Site to retrieve the data
dirProject = 'https://osf.io';

urlStr = [];
fnames = [];

% Define the URL strings for the individual files

% Raw data from NYU subjects (prior to denoising)
rawNYU =        {'ewtz9', 'gfhhk',  ...   conditions & data s01
                'jzx8u', '5yw3e',   ...   conditions & data s02
                'd8me2', 'zyrr2',   ...   conditions & data s03
                '92jhd', '5p8hh',   ...   conditions & data s04
                's7avt', 'hsbm2',   ...   conditions & data s05
                '7a5nq', 'nsqex',   ...   conditions & data s06
                'krdff', 'eycyu',   ...   conditions & data s07
                '6cbqr', 'tjj3v'};      % conditions & data s08

% Denoised data from NYU subjects with 10 PCs removed
denoised10NYU = {'j6tzv', 'md8uu',  ...  denoised BB & SL s01
                'ffgev', 'fk8kb',   ...  denoised BB & SL s02
                'uu5b9', 'gw6up',   ...  denoised BB & SL s03
                'tex2u', '8aaxj',   ...  denoised BB & SL s04
                'ndw48', 'esvau',   ...  denoised BB & SL s05
                'xpc78', 'zh9sr',   ...  denoised BB & SL s06
                '8mrb5', 'sgeam',   ...  denoised BB & SL s07
                'uspmn', 'q8tpt'};     % denoised BB & SL s08

% Denoised data from NYU subjects with all possible PCs removed (1-75)
denoisedAllNYU = {'mf6nx','cw4fr',  ...   denoised BB & SL s01
                '3suq7', 'u77bm',   ...   denoised BB & SL s02
                '9hkns', '74vy5',   ...   denoised BB & SL s03
                '8m4vq', 'h53z8',   ...   denoised BB & SL s04
                's6m42', 'auzuh',   ...   denoised BB & SL s05
                'f4r6m', 'ach36',   ...   denoised BB & SL s06
                'hveh9', 'cpf9x',   ...   denoised BB & SL s07
                'gvpky', '9es6q'};      % denoised BB & SL s08

% Denoised data from NYU subjects with CALM as preprocessing step (see
% Adachi, Shimogawara, Higuchi, Haruta, Ochiai (2001) IEEE T Appl Supercon)
CALMNYU =       {'584hr', 'd92qf', ...    conditions & data s21
                'buqgg', '226gd', ...     conditions & data s22
                'g8cjm', 'xnmqq', ...     conditions & data s23
                '9y9ed', 'd92wr', ...     conditions & data s24
                '554gw', 'vnzjy', ...     conditions & data s25
                'kncfc', 'bqgz2', ...     conditions & data s26
                'yzegk', 'cp24s', ...     conditions & data s27
                'xtu4k', '8hgpy'};      % conditions & data s28

% Denoised data from NYU subjects with TSPCA as preprocessing step (see de
% Cheveigne & Simon (2007) JNeuroci Methods)
TSPCANYU =      {'6fwkd', 'wghxr',  ...   conditions & data s29
                '9jpyr', 'tj3cc',   ...   conditions & data s30
                'uw3ap', '2pexw',   ...   conditions & data s31
                'r5re8', 'qj4d2',   ...   conditions & data s32
                'th6xy', '8g3m4',   ...   conditions & data s33
                '8r97t', '6c947',   ...   conditions & data s34
                '6ufrt', 'y28mr',   ...   conditions & data s35
                'z5kq7', 'b3fs8'};      % conditions & data s36

% Denoised data from NYU subjects using 6 types of control analyses
controlNYU =    {'43p8e', 's7qam', 'suj6e', 'avqwf', 's2xsa', 'ugj6g',   ... Control 1-6 s01
                '6bfh6', '8xk5u', 'y4k2r', 'm77zu', '595cu', '364bq',    ... Control 1-6 s02
                'pevc7', '4r9xy','d827m', '5rw4y', 'kdfu8', 'gxgja',     ... Control 1-6 s03
                '3ctyg', '566es', '2f74r', '2wgs9', 'nfhea', 'ezr3e',    ... Control 1-6 s04
                'uz55w', 'accqa','msxjr','b39gy', 'n7nyg', 'zwgeb',      ... Control 1-6 s05
                'bj495', '4d5kr', 'cv4fp', 'vhtfy', '3nawg', 'wk38v',    ... Control 1-6 s06
                'pyvw5', 'swqhg', '4h8um', 'hfxz9', '8ykes', 'aj2hy',    ... Control 1-6 s07
                '9p9hy', 'h2czn', 'nn64m', '7udf5', '5jraa', 'wx3de'};     % Control 1-6 s08

% Raw data from CiNet subjects (prior to denoising)
rawCiNET =      {'k3vub','n7j57',   ... % conditions & data s09
                'xdt9t', 'xz9d2',   ... % conditions & data s10
                'tj6wn', '9znqn',   ... % conditions & data s11
                'evkmk', 'suu3t'};      % conditions & data s12

% Denoised data from CiNet S9-12 subjects using TSSS as preprocessing step
% are defined as new subject numbers S14, S16, s18, S20. (S13, S15, S17, S19) are
% subject S9, S10, S11, S12 preprocessed with SSS algorithm (Taulu & Simola,
% 2006, Phys Med Biol), this algorithm is highly similar to tSSS (Taulu &
% Hari, 2009, HBM). In the paper we only use the more recent tSSS algorithm
% as a preprocessing step.
% 
TSSSCiNET =     {'7xj4j', '5vgh7',  ... conditions & data s14
                '', '',             ... left out SSS data (s15)
                'csk7e', 'u9xsy',   ... conditions & data s16
                '', '',             ... left out SSS data (s17)
                'w7d7n', 'fq4t2',   ... conditions & data s18
                '', '',             ... left out SSS data (s19)
                'fwvcj', 'ctc99'};      % conditions & data s20

% Denoised data from NYU subjects sweeped across the 3D space parameter space
% varying PCs projected out, varying size of noise pool, varying amount of
% epochs denoised at a time
SF1 =           {'bu9nx', ... SF1 data s01
                'g4aqn', ... SF1 data s02
                'x6nzm', ... SF1 data s03
                'kre9n', ... SF1 data s04
                'fsjwr', ... SF1 data s05
                'tg8xq', ... SF1 data s06
                'bc9hx', ... SF1 data s07
                'tzshp'};  % SF1 data s08

SF2 =           {'8q6tg', ... SF2 data s01
                'vxy2t', ... SF2 data s02
                'm9zqp', ... SF2 data s03
                '2wnp9', ... SF2 data s04
                'dqtu6', ... SF2 data s05
                'bjvd2', ... SF2 data s06
                '6wzkp', ... SF2 data s07
                'b96pd'};  % SF2 data s08

SF3 =           {'79vfs', 'sv2yp'};  % SF3 noisedata & noisepcs s01

% Synthetic data set with BB (broadband) or without broadband
SF4_withBB =    {'4d9jg', 'bf58n', ... conditions & data s99
                'tcx4d', 'a29ve', ... denoised BB & SL s99
                '8mh5c', 'd4uhx', ... denoised Full BB & SL s99
                'jfr4z'};           % denoised timeseries s99

SF4_withoutBB = {'28epj', 'qgesd', ... conditions & data s99
                'nes4t', 'wj37p', ... denoised BB & SL s99
                '39zmb', '5u9aj', ... denoised Full BB & SL s99
                '38jsd'};           % denoised timeseries s99

% Concatenate the url's where we continu to count subject numbers
raw = cat(2,rawNYU,rawCiNET,TSSSCiNET,CALMNYU,TSPCANYU);

% Get postfix of data file
for ii = 1:length(whichDataTypes)
    switch lower(whichDataTypes{ii})
        case 'raw'
            urlStr = [urlStr raw]; %#ok<*AGROW>
            fnames = [fnames, {'_conditions'},{'_sensorData'}];
        case 'denoised 10 pcs'
            urlStr = [urlStr denoised10NYU];
            fnames = [fnames, {'_denoisedData_bb'},{'_denoisedData_sl'}];
        case 'denoised all pcs'
            urlStr = [urlStr denoisedAllNYU];
            fnames = [fnames, {'_denoisedData_full_bb'},{'_denoisedData_full_sl'}];
        case 'controls'
            urlStr = [urlStr controlNYU];
            fnames = [fnames, {'_denoisedData_control'}];
        case 'sf1'
            urlStr = [urlStr SF1];
            fnames = [fnames, {'_denoisedData_NCPSvsNoisePool_bb'}];
        case 'sf2'
            urlStr = [urlStr SF2];
            fnames = [fnames, {'_denoisedData_varyEpochLength_NrPCs_bb'}];
        case 'sf3'
            urlStr = [urlStr SF3];
            fnames = [fnames, {'_noisedata'},{'_noisepcs'}];
        case 'sf4'
            if ~any(whichSubjects==99); error('Can only download SF4 data for synthetic data (subject 99)'); end
            urlStr = [urlStr SF4_withBB SF4_withoutBB];
            fnames = [fnames, repmat([{'conditions'},{'_sensorData'}, ...
                {'_denoisedData_bb'},{'_denoisedData_sl'}, ...
                {'_denoisedData_full_bb'},{'_denoisedData_full_sl'}, ...
                {'_denoisedts'}],1,2)];
            % Make folder in case it doesn't exist
            if ~exist(fullfile(nppRootPath, 'analysis', 'data','s99_withoutBB'),7)
                    mkdir(fullfile(nppRootPath, 'analysis', 'data','s99_withoutBB')); end
    end
end

% Read / write the sample data
for s = whichSubjects
    
    fprintf('Downloading subject %d .\n',s);
    
    % If we deal with synthetic dataset
    if s == 99
        for f = 1:length(fnames)
            fname = sprintf('s%02d%s.mat', s, fnames{f});
            
            % Create path to write file
            writePth = fullfile(savePth, fname);
            
            if any(strcmp(urlStr{f},SF4_withoutBB)) 
                savePth = fullfile(nppRootPath, 'analysis', 'data','s99_withoutBB');
            end
            
            readPth  = fullfile(dirProject, urlStr{f}, '?action=download&version=1');
            websave(writePth,readPth);    
        end
        
        
    else  % All other subjects
        
        for f = 1:length(fnames)
            if ~strcmp(fnames{f},'_denoisedData_control')
                fname = sprintf('s%02d%s.mat', s, fnames{f});
                
                % Create path to write file
                writePth = fullfile(savePth, fname);
                
                % Try to download the second version, if exists, else download
                % version 1.
                try
                    readPth  = fullfile(dirProject, urlStr{f}, '?action=download&version=2');
                    websave(writePth,readPth);
                catch
                    warning('(nppDownloadsampledata): Cannot download version 2 of file, so will default to version 1')
                    readPth  = fullfile(dirProject, urlStr{f}, '?action=download&version=1');
                    websave(writePth,readPth);
                end
                
            else
                for ncontrol = 1:6
                    fname = sprintf('s%02d%s%d_bb.mat', s, fnames{f}, ncontrol);
                    readPth  = fullfile(dirProject, urlStr{ncontrol+(6*(s-1))}, '?action=download&version=2');
                    writePth = fullfile(savePth, fname);
                    websave(writePth,readPth);
                end
            end
        end
    end
end




fprintf('Downloading is done!\n');

return




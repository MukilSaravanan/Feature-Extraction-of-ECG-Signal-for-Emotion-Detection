 classdef Feature_extracter < handle
    
    properties (SetAccess = private)       
        %General
        name 
        BPMglobal {isnumeric}
        BPMlocal {isnumeric}
        len {isnumeric}
        Fs {isnumeric}
        signal {ismatrix}        
        %Numerics
        SDNN {isnumeric}
        rmsRR {isnumeric}
        NN50 {isnumeric}
        pNN50 {isnumeric}
        IBImean {isnumeric}
        IBIrange {ismatrix}
        %Poincare Map
        SD1 {isnumeric}
        SD2 {isnumeric}
        S_area {isnumeric}
        %Power Values
        lfPower {isnumeric}
        hfPower {isnumeric}
        lowFreqHighFreqRatio {isnumeric}
        TotPower {isnumeric}
        LFperc {isnumeric}
        HFperc {isnumeric}
        %Segments
        SegmentCount {isnumeric}
        SegmentLength {isnumeric}
        Segments {ismatrix}
    end
    
    properties (Access = private)
        %General
        timeVal {ismatrix}
        pks {ismatrix}
        locs {ismatrix}
        %HRV
        rr {ismatrix}
        pkTimes {ismatrix}
        rrA {ismatrix} %Backup
        pkTimesA {ismatrix} %Backup
        reject {ismatrix}
        HRV_R {ismatrix} %Resampled RR
        tHRV {ismatrix}
        LHRV {isnumeric}
        FsHRV {isnumeric}
        %Freq Analysis
        fftHRV {ismatrix}
        p1 {ismatrix}
        %Freq Bands Analsysis
        hf {ismatrix}
        p1HF {ismatrix}
        lf {ismatrix}
        p1LF {ismatrix}
        %Segments
        segmentBoundaries {ismatrix}
    end
    
    methods 
        %% General Calculation
        function obj = Feature_extracter(ECGSignal,SamplingFrequency,Name)
            %Constructor, (Supply Signal and Initial Sampling Frequency)
            %Creates ECG class with many methods for analysis
            
            if nargin >= 2
                obj.signal=ECGSignal;
                obj.Fs=SamplingFrequency;
                
                obj.len = length(obj.signal);
                T = 1/obj.Fs;
                obj.timeVal = (0:obj.len-1)'*T;
                
           
                obj.name='ECG';
            end
            
            if nargin == 3
                obj.name = Name;
            end

        end
        
        function init(obj,varargin)
            % Eliminates Offset, Detrends, Detects Peaks, Gets RR data and
            % filters it
            % Paramters are:
            % Disabling sections: 'offSetElim', 'deTrend'
            % Rejection: 'rejectPlot', 'STDscalar'
            % Peak Detect: 'MinPeakProminence','MinPeakHeight','MinPeakDistance'
            % Detrending: 'ORD', 'FL'
            
            p = inputParser;
            p.KeepUnmatched = true;
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            addParameter(p,'offSetElim','on',checkPlot)
            addParameter(p,'deTrend','on',checkPlot)
            parse(p,varargin{:});
            
            if strcmpi(p.Results.offSetElim,'on')
                obj.offsetEliminate
            end
            if strcmpi(p.Results.deTrend,'on')
                obj.deTrending(varargin{:})
            end
            obj.peakDetect(varargin{:})
            obj.variabilityReject(varargin{:})
            obj.calculateHRVNumerics
            obj.resampleHRV
            obj.freqAnalysisHRV
        end
        
        function resample(obj,FsR,varargin)
            % Resamples the ECG at a given resampling frequency
            % Parameters are Plot (on/off), ylim, xlim
            TR = 1/FsR;
            p = inputParser;
            validLim = @(x) (ismatrix(x) && length(x) == 2);
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            addParameter(p,'ResamplePlot','off',checkPlot)
            addParameter(p,'xlim',[0, max(obj.timeVal)],validLim)
            addParameter(p,'ylim',[min(obj.signal)*1.2, max(obj.signal)],validLim)
            parse(p,varargin{:});
            
            if (strcmp(p.Results.ResamplePlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','m')
                hold on
            end
            
            obj.signal = resample(obj.signal,FsR,obj.Fs,FsR);
            obj.len = length(obj.signal);
            
            obj.timeVal = (0:obj.len-1)*TR;
            
            if (strcmp(p.Results.ResamplePlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','k')
                legend(['Original ECG, sampled at ',num2str(obj.Fs),' Hz'],...
                       ['Resampled ECG at ',num2str(FsR),' Hz'])
                hold off
                
            end
            
            obj.Fs = FsR;
        end
        
        function trim(obj,range,varargin)
            %Trims ECG signal outside a given range (seconds)
            %Params are TrimPlot, trimColour, sigColour, xlim, ylim
            
            p = inputParser;
            validLim = @(x) (ismatrix(x) && length(x) == 2);
            checkPlot = @(x) (strcmpi(x,'on') || strcmpi(x,'off'));
            addParameter(p,'TrimPlot','off',checkPlot)
            addParameter(p,'trimColour',[0.2471, 0.2627, 0.2784, 0.15])
            addParameter(p,'sigColour',[0, 0.3569, 0.6588])
            addParameter(p,'xlim',[0, max(obj.timeVal)],validLim)
            addParameter(p,'ylim',[min(obj.signal)*1.2, max(obj.signal)],validLim)
            parse(p,varargin{:});
            
            
            
            cuttings = find(obj.timeVal > max(range) | obj.timeVal < min(range));
            
            if strcmpi(p.Results.TrimPlot,'on')
                obj.plot('Colour',p.Results.trimColour)
                hold on
            end
            maxT = max(obj.timeVal);
            
            
            obj.signal(cuttings) = [];
            obj.timeVal(cuttings) = [];
            obj.len = length(obj.signal);
            
            if strcmpi(p.Results.TrimPlot,'on')               
                obj.plot('Colour',p.Results.sigColour,'xlim',[0 maxT])
                hold off
            end
            
        end
        
        function resetTimeVal (obj)
            %Resets time values to start at 0s, following trimming or
            %segmentation
            offset = min(obj.timeVal);
            obj.timeVal = obj.timeVal - offset;
            
        end
        
        function offsetEliminate(obj)
            %Eliminates any offsets in ECG signal
            offset = mean(obj.signal);
            obj.signal = obj.signal - offset;
        end
        
        function deTrending(obj,varargin)
            %Removes trend
            %input ORD (order) and FL (Float length)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'ORD',3)
            addParameter(p,'FL',9999)
            validLim = @(x) (ismatrix(x) && length(x) == 2);
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            addParameter(p,'trendPlot','off',checkPlot)
            addParameter(p,'xlim',[0, max(obj.timeVal)],validLim)
            addParameter(p,'ylim',[min(obj.signal), max(obj.signal)],validLim)
            
            parse(p,varargin{:});
            
            if (strcmp(p.Results.trendPlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','m')
                hold on
            end
            
            trend = sgolayfilt(obj.signal, p.Results.ORD, p.Results.FL);
            obj.signal = obj.signal-trend;
            
            if (strcmp(p.Results.trendPlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','b')                
                plot(obj.timeVal,trend,'k')
                xlim(p.Results.xlim)
                ylim(p.Results.ylim)
                hold off
            end
        end
        
        function filter(obj,SecondOrderSection,Gain,varargin)
            % Filters the signal using a supplied SecondOrderSection and
            % Gain
            % Parameters are: FiltPlot,fftPlot, xlim and ylim
            p = inputParser;
            validLim = @(x) (ismatrix(x) && length(x) == 2);
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            addParameter(p,'filtPlot','off',checkPlot)
            addParameter(p,'fftPlot','off',checkPlot)
            addParameter(p,'xlim',[0, max(obj.timeVal)],validLim)
            addParameter(p,'ylim',[min(obj.signal)*1.2, max(obj.signal)],validLim)
            parse(p,varargin{:});
            
            signalTemp = obj.signal;
            
            if (strcmp(p.Results.fftPlot,'on'))
                obj.plotSignalFFT('Colour','m','ylim',p.Results.ylim)
                hold on
            end
            
            if (strcmp(p.Results.filtPlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','m')
                hold on
            end
            
            obj.signal = filtfilt(SecondOrderSection,Gain,obj.signal);
            
            if (strcmp(p.Results.filtPlot,'on'))
                obj.plot('xlim',p.Results.xlim,'ylim',p.Results.ylim,'Colour','k')
                legend('Unfiltered ECG','Filtered ECG')
                hold off
            end  
            
            if (strcmp(p.Results.fftPlot,'on'))
                obj.plotSignalFFT('Colour','k','ylim',p.Results.ylim)
                legend('Unfiltered Signal','Filtered Signal')
                hold off
            end
            
        end
        
        function peakDetect(obj,varargin)
            % Finds Peaks of ECG, enter parameter pair values of
            % MinPeakProminence, MinPeakDistance, MinPeakHeight
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'MinPeakProminence',0.1)
            addParameter(p,'MinPeakDistance',0.4)
            addParameter(p,'MinPeakHeight',0.5)
            
            validLim = @(x) (ismatrix(x) && length(x) == 2);
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            addParameter(p,'peakDetectPlot','off',checkPlot)
            addParameter(p,'xlim',[0, max(obj.timeVal)],validLim)
            addParameter(p,'ylim',[min(obj.signal), max(obj.signal)],validLim)
            
            parse(p,varargin{:});
            
            PP=p.Results.MinPeakProminence;
            PD=p.Results.MinPeakDistance;
            PH=p.Results.MinPeakHeight;
            
            [obj.pks,obj.locs] = findpeaks(obj.signal,obj.timeVal,'MinPeakProminence',PP,...
                'MinPeakDistance',PD,'MinPeakHeight',PH);
            
            if (strcmp(p.Results.peakDetectPlot,'on'))
            plot(obj.timeVal,obj.signal,'k',obj.locs,obj.pks,'rd')
            xlim(p.Results.xlim)
            ylim(p.Results.ylim)
            grid on
            end
            
            
            obj.BPMglobal = length(obj.locs)*100*60/obj.len;
            obj.rr = diff(obj.locs);
            obj.pkTimes = obj.locs(1:end-1);
            obj.BPMlocal = 60/mean(obj.rr);
        end
        
        function variabilityReject(obj,varargin)
            % Rejects all noise from RR-Interval data
            % Input args are 'rejectplot','on'/'off' or 'STDscalar',(value)
            
            p = inputParser;
            p.KeepUnmatched = true;
            checkPlot = @(x) (strcmp(x,'on') || strcmp(x,'off'));
            checkScalar = @(x) (isnumeric(x) && x > 0 && x < 4);
            addParameter(p,'RejectPlot','off',checkPlot)
            addParameter(p,'STDscalar',1.96,checkScalar)
            addParameter(p,'Reject','on',checkPlot)
            parse(p,varargin{:});
            
            if (strcmp(p.Results.Reject,'on'))
                
                pkMean = mean(obj.rr);
                pk2STD = std(obj.rr) * p.Results.STDscalar;
                rejectMax = pkMean + pk2STD;
                rejectMin = pkMean - pk2STD;
                
                ArrReject = find((obj.rr > rejectMax) | (obj.rr < rejectMin));
                
                if (strcmp(p.Results.RejectPlot,'on'))
                    
                    obj.reject = zeros(length(ArrReject),2);
                    
                    for i = 1:length(ArrReject)
                        obj.reject(i,1) = obj.rr(ArrReject(i));
                        obj.reject(i,2) = obj.pkTimes(ArrReject(i));
                    end
                    
                    plot(obj.pkTimes,obj.rr,'b'), hold on
                    scatter(obj.reject(:,2),obj.reject(:,1),'ro')
                    yline(rejectMax,'--r')
                    yline(rejectMin,'--r')
                    ylabel('RR Intervals(s)')
                    xlabel('Time (s)')
                    title('RR Interval Series Showing Rejected Points')
                    hold off
                    
                end
                
                obj.rrA=obj.rr;
                obj.pkTimesA=obj.pkTimes;
                obj.rr(ArrReject) = [];
                obj.pkTimes(ArrReject) = [];
                
                obj.BPMlocal = 60/mean(obj.rr);
                
            end
        end
        
        function resetRR(obj)
            obj.rr=obj.rrA;
            obj.pkTimes=obj.pkTimesA;
            obj.BPMglobal = mean(obj.rr)*60;
        end
        
        function calculateHRVNumerics(obj)
            % Calculates SDNN, RMSSD, NN50, and pNN50 values
            differences = diff(obj.rr);
            obj.SDNN = std(obj.rr*1000);
            obj.rmsRR = sqrt(mean(differences.^2))*1000;
            
            obj.NN50 = 0;
            
            for i = 1:length(differences)
                if (abs(differences(i)*1000) > 50)
                    obj.NN50 = obj.NN50 + 1;
                end
            end
            
            obj.pNN50 = obj.NN50*100 / length(differences);
            
            % Poincare Map
            
            pmX = obj.rr*1000;
            pmX(end)=[];
            
            pmY = obj.rr*1000;
            pmY(1)=[];
            
            obj.SD1 = std(pmX - pmY);
            obj.SD2 = std(pmX + pmY);
            
            % Beat Interval
            obj.IBImean = mean(obj.rr*1000);
            obj.IBIrange = [max(obj.rr*1000) min(obj.rr*1000)];
        end
        
        function resampleHRV(obj,FsHRV)
            if nargin < 2 || isempty(FsHRV)
                obj.FsHRV = 4;
            else
                obj.FsHRV = FsHRV;
            end
            
            THRV = 1/obj.FsHRV;
            obj.HRV_R = resample(obj.rr,obj.pkTimes,obj.FsHRV);
            obj.LHRV = length(obj.HRV_R);
            obj.tHRV = (0:obj.LHRV-1)'*THRV;
        end
        
        function freqAnalysisHRV(obj)
            % Analyses High, Low and Total Freqeuncy and Power of RR
            % interval series
            warning('off','MATLAB:colon:nonIntegerIndex')
            yF = fft(obj.HRV_R);
            p2 = abs(yF/obj.LHRV);
            obj.p1 = p2(1 : obj.LHRV/2 + 1);
            obj.p1(2 : end-1) = 2 * obj.p1(2 : end-1);
            
            obj.fftHRV = obj.FsHRV * (0 : (obj.LHRV/2)) / obj.LHRV;
            
            % LF
            
            NOTlf = find(obj.fftHRV < 0.04 | obj.fftHRV > 0.15);
            obj.lf = obj.fftHRV;
            obj.lf(NOTlf) = [];
            
            obj.p1LF = obj.p1;
            obj.p1LF(NOTlf) = [];
            
            obj.lfPower = sum(obj.p1LF)*1000;
            % HF
            
            NOThf = find((obj.fftHRV > 0.4) | obj.fftHRV < max(obj.lf)); %| obj.fftHRV < 0.15);
            obj.hf = obj.fftHRV;
            obj.hf(NOThf) = [];
            
            obj.p1HF = obj.p1;
            obj.p1HF(NOThf) = [];
            
            obj.hfPower = sum(obj.p1HF)*1000;
            
            obj.lowFreqHighFreqRatio = obj.lfPower/obj.hfPower;
            obj.TotPower = obj.lfPower + obj.hfPower;
            obj.LFperc = obj.lfPower*100/(obj.lfPower + obj.hfPower);
            obj.HFperc = 100-obj.LFperc;
            warning('on','MATLAB:colon:nonIntegerIndex')
        end
        
        function rename(obj,newName)
            
            obj.name = newName;
        
        end 
        
        function features=get_features(obj)
            bpmglobal=obj.BPMglobal;
            bpmlocal=obj.BPMlocal;
            sdnn=obj.SDNN;
            rmsrr=obj.rmsRR;
            nn50=obj.NN50;
            pnn50=obj.pNN50;
            ibimean=obj.IBImean;
            ibirange=obj.IBIrange;
            sd1=obj.SD1;
            sd2=obj.SD2;
            lfpower=obj.lfPower;
            hfpower=obj.hfPower;
            lfhfratio=obj.lowFreqHighFreqRatio;
            totpower=obj.TotPower;
            lfperc=obj.LFperc;
            hfperc=obj.HFperc;
            features=[bpmglobal,bpmlocal,sdnn,rmsrr,nn50,pnn50,ibimean,ibirange,sd1,sd2,lfpower,hfpower,lfhfratio,totpower,lfperc,hfperc];
            
        end
    end
    
    methods (Access = private)
        %% Private Methods
        
        function incrementTimeVal(obj,increment)
            
            obj.timeVal = obj.timeVal + increment;
            
        end  
        
    end
end

clear ;
clc;

%---------------------The Synthetic Signal ---------------------------------
N = 120;                 % Samples per cycle
f0 = 60;                 % Fundamental Frequency of the signal in Hz
fs = N * f0;             % Sampling frequency
T = 1/fs;                % Sample time or rate delta T
t = (0:T:(1/f0)*2);           % Time vector
amplitude_drift = 10/100;


h = 5;                  %no of harmonics
AmpPhaseValues    =zeros(1,2*h);  %Amplitude and Phase for each component

PureSignal      = zeros(length(t),1);
Measured_Y      = zeros(length(t),1);

FundPhase = 5*pi/180;           %5 degrees
FundAmp   = 1;                  % 1 p.u

FundAmpVector = ones(length(t),1);
%  60 Hz sinusoid and 4 harmonics
for i=1:h
    AmpPhaseValues(2*i-1)   = (FundAmp/i);
    AmpPhaseValues(2*i) = i*FundPhase;
    PureSignal =PureSignal + (AmpPhaseValues(2*i-1)*sin(2*pi*i*f0*t +AmpPhaseValues(2*i)))';
end
%Infected with noise
signalTonoisedB = 120;
Measured_Y =awgn(PureSignal,signalTonoisedB,'measured');

NoOfWeights = 2*h;      % total weights of all component
%----------Using the Case B estimation model and harmonics---------
InputSignalVector_X = zeros(length(t),NoOfWeights);
for n =1:length(t)
    for i =1:h
        InputSignalVector_X(n,2*i-1)= sin(2*pi*i*f0*t(n));
        InputSignalVector_X(n,2*i)= cos(2*pi*i*f0*t(n));
    end
end
%--------------------------------------------------------------------------

%--------------------LMS Conj Algorithm and Computations of Phasors--------
noOfSamples = (NoOfWeights)*3;

Algo = BlockLMSConjAlgorithm(noOfSamples,1,InputSignalVector_X,Measured_Y);
Algo.Process(InputSignalVector_X,Measured_Y);
EstimatedFundAmplitude    = zeros(length(t),1);
TVEofHarmonics             = zeros(length(t),NoOfWeights/2);

for n =1:length(t)
    if((n+noOfSamples-1)<length(t))
        weights = Algo.EvolvedWeightVectors(n,:);
        for i =1:h
            estAmp   = sqrt(weights(2*i-1)^2+weights(2*i)^2);     %Amplitude of ith harmonic
            estPhase =  atan(weights(2*i)/weights(2*i-1));       %Phase of ith harmonic
            percErrorAmp = 100*(estAmp-AmpPhaseValues(2*i-1))/AmpPhaseValues(2*i-1);
            phaseErrorinDegrees = (180/pi)*((AmpPhaseValues(2*i)-estPhase)/.573);
            TVEofHarmonics(n,i) = 0.01*sqrt(percErrorAmp^2+phaseErrorinDegrees^2);
            
            if(i== 1)
                 EstimatedFundAmplitude(n+noOfSamples-1) = estAmp;
            end
        end
    end
    
end
%-------------------------------------------------------------------------
% Plot fundamental amplitude

hold on;
% plot(t,FundAmpVector,t,EstimatedFundAmplitude,GetRandomLineColor_StylesAndMarker());
% title('Estimated Fundamental Amplitude');
% ylabel('A(p.u)');
% xlabel('t(s)');
% legend('true value','estimated');

for i =1:1
    plot(t,TVEofHarmonics(:,i),GetRandomLineColor_StylesAndMarker())
    title('Static Result tests')
    s{i} = strcat('Harmonic',' ',int2str(i));
    ylabel('TVE(%)')
    xlabel('t(s)')
end
legend(s{:});
hold off;

(.25*N/fs)
max(TVEofHarmonics(:,1))



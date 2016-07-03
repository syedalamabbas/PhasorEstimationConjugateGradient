clear ;
clc;

%---------------------The Synthetic Signal ---------------------------------
N = 120;                      % Samples per cycle
f0 = 60;                 % Fundamental Frequency of the signal in Hz
fs = N * f0;             % Sampling frequency
T = 1/fs;                    % Sample time or rate delta T
t = (0:T:.08);               % Time vector
phase_drift = pi/18;

FundAmpl        = zeros(1,length(t));
FundPhase       = zeros(1,length(t));
FundFrequency   = zeros(1,length(t));
Measured_Y      = zeros(length(t),1);
PureSignal      = zeros(length(t),1);
%  60 Hz sinusoid and a second weaker harmonic
for n =1:length(t)
    FundFrequency(n) = f0;
    FundAmpl(n)=1;%p.u
    if((0.02 <= t(n))&& (t(n) <= 0.06))
        FundPhase(n)= (0)+phase_drift;
    else
        FundPhase(n)= (0);
    end
    PureSignal(n)=FundAmpl(n)*sin(2*pi*FundFrequency(n)*t(n)+ FundPhase(n));
end

FundPhase = (180/pi)*FundPhase;

Measured_Y = PureSignal;
%----------Using the Case B estimation model and using drift =1 second harmonic only---------
NoOfWeights = 2;
InputSignalVector_X = zeros(length(t),NoOfWeights);
for n =1:length(t)
    InputSignalVector_X(n,1)= sin(2*pi*f0*t(n));
    InputSignalVector_X(n,2)= cos(2*pi*f0*t(n));
end
%--------------------------------------------------------------------------

%--------------------LMS Conj Algorithm and Computations of Phasors--------
noOfSamples = NoOfWeights*6;
EstimatedFundPhase    = zeros(1,length(t));
Algo = BlockLMSConjAlgorithm(noOfSamples,NoOfWeights,InputSignalVector_X,Measured_Y);
Algo.Process(InputSignalVector_X,Measured_Y);
%--------------------------------------------------------------------------

%---------------------Computations of Phasors----------------------------
for n =1:length(t)
    if((n+noOfSamples-1)<length(t))
        weights = Algo.EvolvedWeightVectors(n,:);
        EstimatedFundPhase(n+noOfSamples-1) = (180/pi)*atan(weights(2)/weights(1));
    end
end
%-------------------------------------------------------------------------

% Plot fundamental amplitude
figure,plot(t,FundPhase,t,EstimatedFundPhase,GetRandomLineColor_StylesAndMarker())
title('Dynamic response for the phase-angle step')
axis([0.0 0.08 -6 16 ])
ylabel('Phase(degrees)')
xlabel('t(s)')
legend('true value','estimated','estimated LES');



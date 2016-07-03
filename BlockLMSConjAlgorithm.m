classdef BlockLMSConjAlgorithm < handle
    %BlockLMSCONJALGORITHM determines the successful optimal weights for a given
    %linear filter using conjugate search directions in block update
    
    properties
        block_size;          % represents the block size of Input to be used
        length_L;            % represents the length of weight vector
        weights_W;           % represents the weights W vector
        MSE;                 % Error vector Generated in each Iteration
        error;               % Gives the Instantaneous error vector at each pattern
        NoOfPatterns;        % Number of Input vector sequences x(1),x(2)...
        directionVector_s;   % The direction vector for conjugate based update
        residualdirection_r; % The residual steepest direction vector for conjugate based update
        Beta;                % scaling factor for linear combination of residual and direction
        Eta;                 % optimally derived scaling parameter for weight update
        EvolvedWeightVectors;% The weight vectors which are evolved over time
        maxIterations;       %The Algorithm exits if when the input signals are exhausted or maxiterations are reached
    end
    
    methods
        function Algo = BlockLMSConjAlgorithm(block_size,iterations,InputSignalVector_X,Measured_Y)
            Algo.block_size = block_size;
            Algo.maxIterations = iterations;
            Algo.Initialize(InputSignalVector_X,Measured_Y);
        end  %Algorithm Constructor
        
        %--------------------------------------------------------------------------
        %-------This is where the main processing of the Algorithm happens---------
        function W = Process(Algo,InputSignalVector_X,Measured_Y)% Calculates the Current weights using this Algorithm
                                                            % STEP 1
            InputBlock = zeros (Algo.block_size,Algo.length_L);
            OutputBlock = zeros (1,Algo.block_size);
            
            window = 1;
            while( (Algo.NoOfPatterns - window*Algo.block_size) >= 0) 
                startIndex = (window-1)*Algo.block_size +1;
                stopIndex = window*Algo.block_size;
                InputBlock(1:Algo.block_size,:) = InputSignalVector_X(startIndex:stopIndex,:);
                OutputBlock(1:Algo.block_size) = Measured_Y(startIndex:stopIndex);
                window = window +1;
                
                iterationsCount = Algo.maxIterations;           % For each block iterate given times
                while(iterationsCount > 0)
                    iterationsCount = iterationsCount-1;
%                     [RandInput_X,RandDesired_D] = RandomizePairPatterns(InputBlock,OutputBlock');
%                     InputBlock = RandInput_X;
%                     OutputBlock = RandDesired_D';
%                     
                    Algo.error(stopIndex,:) = (OutputBlock' - InputBlock*Algo.weights_W')';
                    Algo.MSE(stopIndex)     = (2/Algo.block_size)*Algo.error(stopIndex,:)*Algo.error(stopIndex,:)';    %Normalized Error Energy function
                    newResidualdirection = (2/Algo.block_size)*(Algo.error(stopIndex,:)*InputBlock);                 % STEP 2
                    
                    tempNumerator        = newResidualdirection *(newResidualdirection-Algo.residualdirection_r)';
                    tempDenominator      = Algo.residualdirection_r * Algo.residualdirection_r';
                    if(tempDenominator ~= 0)
                        Algo.Beta            = max(tempNumerator/tempDenominator,0);                     % STEP 3
                    else
                        Algo.Beta            = 0;
                    end
                    Algo.residualdirection_r = newResidualdirection;
                    Algo.directionVector_s   = Algo.residualdirection_r + Algo.Beta * Algo.directionVector_s;% STEP 4
                    
                    tempNumerator = Algo.error(stopIndex,:)*InputBlock*Algo.directionVector_s'+Algo.directionVector_s*InputBlock'*Algo.error(stopIndex,:)';
                    tempDenominator = 2*(Algo.directionVector_s*(InputBlock'*InputBlock)*Algo.directionVector_s');
                    if(tempDenominator ~= 0)
                        Algo.Eta  =  tempNumerator / tempDenominator;                                        % STEP 5
                    end
                    Algo.weights_W  = Algo.weights_W + Algo.Eta*Algo.directionVector_s;                      % STEP 6
                end
                for i=1:Algo.block_size
                    Algo.EvolvedWeightVectors(startIndex+i-1,:)   =  Algo.weights_W;
                end
                W = Algo.weights_W;                   %The output,if used is just the Weight which are so far best estimates
            end
                
%             for n = 1:Algo.NoOfPatterns
%                 iterationsCount = Algo.maxIterations;
%                 while(iterationsCount > 0)
%                     iterationsCount = iterationsCount-1;
%                     if(n-1 <= (Algo.NoOfPatterns-Algo.block_size))
%                         for i = 1:Algo.block_size
%                             InputBlock(i,:) = InputSignalVector_X(n-1+i,:);
%                             OutputBlock(i)  = Measured_Y(n-1+i);
%                         end
%                         
%                         [RandInput_X,RandDesired_D] = RandomizePairPatterns(InputBlock,OutputBlock');
%                         InputBlock = RandInput_X;
%                         OutputBlock = RandDesired_D';
%                         
%                         Algo.error(n,:) = (OutputBlock' - InputBlock*Algo.weights_W')';
%                         Algo.MSE(n)     = (2/Algo.block_size)*Algo.error(n,:)*Algo.error(n,:)';    %Normalized Error Energy function
%                         newResidualdirection = (2/Algo.block_size)*(Algo.error(n,:)*InputBlock);                 % STEP 2
%                         
%                         if(Algo.block_size > 1)
%                             if(n-1 >= 0)
%                                 tempNumerator        = newResidualdirection *(newResidualdirection-Algo.residualdirection_r)';
%                                 tempDenominator      = Algo.residualdirection_r * Algo.residualdirection_r';
%                                 if(tempDenominator ~= 0)
%                                     Algo.Beta            = max(tempNumerator/tempDenominator,0);                     % STEP 3
%                                 else
%                                     Algo.Beta            = 0;
%                                 end
%                             end
%                         end
%                         %                     fprintf( '\nThe Beta at iteration %g',n);
%                         %                     Algo.Beta
%                         Algo.residualdirection_r = newResidualdirection;
%                         Algo.directionVector_s   = Algo.residualdirection_r + Algo.Beta * Algo.directionVector_s;% STEP 4
%                         
%                         tempNumerator = Algo.error(n,:)*InputBlock*Algo.directionVector_s'+Algo.directionVector_s*InputBlock'*Algo.error(n,:)';
%                         tempDenominator = 2*(Algo.directionVector_s*(InputBlock'*InputBlock)*Algo.directionVector_s');
%                         if(tempDenominator ~= 0)
%                             Algo.Eta  =  tempNumerator / tempDenominator;                                        % STEP 5
%                         end
%                         
%                         Algo.weights_W  = Algo.weights_W + Algo.Eta*Algo.directionVector_s;                      % STEP 6
%                         Algo.EvolvedWeightVectors(n,:)   =  Algo.weights_W;
%                     end
%                     %                 fprintf( '\nThe weights at iteration %g',Algo.maxIterations-iterationsCount);
%                     %                 Algo.weights_W
%                 end
%             end
%             
%             W = Algo.weights_W;                   %The output,if used is just the Weight which are so far best estimates
        end
        %--------------------------------------------------------------------------
        
        %-----------------------Initilization--------------------------------------
        function Initialize(Algo,InputSignalVector_X,Measured_Y)
            [pInput,kInput]      = size(InputSignalVector_X) ; % represents the number of Patterns vs Input Dimentionality
            [pOutput,kOutput]    = size(Measured_Y);           % represents the number of Patterns vs Output Dimentionality
            if(pInput ~= pOutput)
                error('\nInput and Output number of patterns mismatch')
            end
            if(kOutput ~= 1)
                error('\nThe output has to be scalar')
            end
            
            Algo.NoOfPatterns           = pInput;
            Algo.length_L               = kInput;
            Algo.Beta                   = 0;
            Algo.Eta                    = 0;
            
            Algo.weights_W              = rand(1,Algo.length_L);
            Algo.directionVector_s      = zeros(1,Algo.length_L);
            Algo.residualdirection_r    = zeros(1,Algo.length_L);
            
            Algo.residualdirection_r    = (Measured_Y(1)* InputSignalVector_X(1,:));
            Algo.directionVector_s      = Algo.residualdirection_r;
            
            Algo.error                  = zeros(Algo.NoOfPatterns,Algo.block_size);
            Algo.MSE                    = zeros(1,Algo.NoOfPatterns);
            
            Algo.EvolvedWeightVectors   = zeros(Algo.NoOfPatterns,Algo.length_L);
            Algo.EvolvedWeightVectors(1,:) =  Algo.weights_W;
        end
        %-------------------------------------------------------------------------
    end %Methods ends
end % Class Ends


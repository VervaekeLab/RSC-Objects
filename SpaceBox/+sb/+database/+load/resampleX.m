function Y = resampleX(X, P, Q)


%     % Resample pupil area and pupil movement
%     P = nOrigSamples(block);
%     Q = nPupilSamples;
   
    POrig = P;

    % Limitation of product of resampling factors.
    if P*Q > 2^31
        if ~mod(P, 2) == 0
            P = (P-1)/2;
        else
            P = P/2;
        end
        
        if ~mod(Q, 2) == 0
            Q = (Q+1)/2;
        else
            Q = Q/2;
        end
    end
    
    Y = resample(X, P, Q);
    
    % Correct for resampling limitation
    if numel(Y) < POrig
        Y(POrig) = 0;
    elseif numel(Y) > POrig
        Y = Y(1:POrig);
    end

end
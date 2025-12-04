function Pq = quantize(Phi_in, qbits)
    % Quantize phases to qbits resolution and magnitudes to a few levels.
    % This is a simple model: quantize phase to 2^qbits levels and keep magnitude as-is.
    %
    % Input:
    %   Phi_in  - NxN diagonal or general matrix (we quantize diagonal entries/angles)
    %   qbits   - integer number of phase quantization bits
    %
    % Output:
    %   Pq - quantized matrix same size as Phi_in
    
    % If Phi_in is diagonal, operate on diag; else operate elementwise angles/mags for diagonal-dominant case.
    [n1,n2] = size(Phi_in);
    Pq = Phi_in;
    levels = 2^qbits;
    % operate on diagonal entries if matrix is diagonal (common for IRS)
    if isequal(Phi_in, diag(diag(Phi_in)))
        ph = angle(diag(Phi_in));
        mag = abs(diag(Phi_in));
        % map phase to [0, 2pi)
        ph = mod(ph, 2*pi);
        ph_idx = round(ph/(2*pi) * (levels-1));
        ph_q = ph_idx/(levels-1) * 2*pi;
        diag_q = mag .* exp(1j*ph_q);
        Pq = diag(diag_q);
    else
        % elementwise fallback: quantize angle only
        mags = abs(Phi_in);
        phs = mod(angle(Phi_in), 2*pi);
        ph_idx = round(phs/(2*pi) * (levels-1));
        ph_q = ph_idx/(levels-1) * 2*pi;
        Pq = mags .* exp(1j*ph_q);
    end
end

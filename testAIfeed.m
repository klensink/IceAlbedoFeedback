clear;
T = [260 290 295 289 295 267];

Tlow = 270;
Thigh = 290;

ai = 0.6;
ao =[0.3200    0.1800    0.1900    0.2350    0.2700    0.3050];
a = zeros(size(ao));
%%
%Classify current temps
i_below = T <= Tlow;
i_middle = T > Tlow & T < Thigh;
i_above = T >= Thigh;

% If below threshold, turn to ice
if nnz(i_below) > 0;
    
    a(i_below) = ai;
end

% If inbetween threshold, calc new albedo
if nnz(i_middle) > 0;
    
    a(i_middle) = ao(i_middle) + (ai - ao(i_middle)).*((T(i_middle)...
                   - Thigh).^2)./((Tlow - Thigh)^2);
               
end

% If above threshold, turn to ice
if nnz(i_above) > 0;
    
    a(i_above) = ao(i_above);
    
end

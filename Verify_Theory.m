clear all;


S  = 100000;
T1 = 10;
T2 = 100;

s = 4.8;
h0 = 50;
w0 = 2.3;

f_thresh = 24;

T_max = max(T1,T2);

for i=1:3
    
    if i==1
        s_h = 2;
        s_w = 0.1;
        d   = 1;
    elseif i==2
        s_h = 0;
        s_w = 0;
        d   = 1;
    elseif i==3
        s_h = 2;
        s_w = 0.1;
        d   = 6;
    end
    
    
    x = zeros(S, T_max, d);
    c = zeros(S, T_max, d);
    h = zeros(S, T_max);
    w = zeros(S, T_max);
    
    c(:,1,:) = 0;
    h(:,1)   = h0;
    w(:,1)   = w0;
    
    for t=1:T_max-1
        if i==1 || i==2
            c(:,t+1,:) = c(:,t,:) + s*sign(randn(S,d));
        elseif i==3
            c(:,t+1,:) = squeeze(c(:,t,:)) + s*randn(S,d);
        end
        h(:,t+1)   = h(:,t) + s_h*randn(S,1);
        w(:,t+1)   = w(:,t) + s_w*randn(S,1);
    end
    
    f = h - w.*vecnorm(x-c,2,3);
    
    if i==1
        
        metric1 = mean(f);
        metric1 = cumsum(metric1(:,1:T1)) ./ (1:T1);
        
        aux = 0;
        metric2 = zeros(1,T1);
        metric2(1) = h0;
        for t=1:T1-1
            
            aux = aux + s/2^(t-1)*(floor(t/2)+1)*nchoosek(t, floor(t/2)+1);
            
            metric2(t+1) = h0 - w0/(t+1)*aux;
            
        end
        
    elseif i==2
        
        metric1 = ceil((h0-f_thresh)/(s*w0))^2;
        
        [~,metric2] = min(f > f_thresh, [], 2);
        metric2     = metric2 - 1;
        
        ii          = sum(f > f_thresh,2) == T2;
        metric2(ii) = T2;
        
        metric2 = mean(metric2);
        
    elseif i==3
        
        metric1 = mean(f);
        metric1 = cumsum(metric1(:,1:T1)) ./ (1:T1);
        
        aux = 0;
        metric2 = zeros(1,T1);
        metric2(1) = h0;
        for t=1:T1-1
            
            aux = aux + sqrt(2)*s*sqrt(t)*gamma((d+1)/2)/gamma(d/2);
            
            metric2(t+1) = h0 - w0/(t+1)*aux;
            
        end
        
    end
    
    [metric1; metric2]
    
end






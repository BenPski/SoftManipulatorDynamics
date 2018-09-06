function out = interpNet(net,v,t,train_interval,des_interval)
    %interpolates the changes in state if this is a change in state
    %predictor
    out = (net([v;t])-t)*des_interval/train_interval+t;
end
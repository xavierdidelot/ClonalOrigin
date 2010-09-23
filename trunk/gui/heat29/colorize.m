function vec=colorize(gain,loss)
    vec=[1 1 1];
    if gain>0&&gain<0.5,vec=[1,1,1-(gain*2)^2];end
    if gain>=0.5,vec=[1,1-((gain-0.5)*2)^2,0];end
    if loss>0&&loss<0.5,vec=[1-(loss*2)^2,1,1];end
    if loss>=0.5,vec=[0,1-((loss-0.5)*2)^2,1];end
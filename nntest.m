load hourno0.mat
T = hourno0; 
net = narnet(1:2,10); 
[Xs,Xi,Ai,Ts] = preparets(net,{},{},T); 
net = train(net,X,T,Xi,Ai); 
view(net) 
Y = net(Xs,Xi,Ai) 
plotresponse(T,Y)
pd = makedist('Normal');

sigma = [1.5, 3,4, 4.5 ,5,6];

cpk = sigma./3;

maxNC = (2*cdf(pd, - (sigma))) * 10^6
minNC = (cdf(pd, - sigma))*10^6;


result = [cpk', sigma', maxNC', minNC'];
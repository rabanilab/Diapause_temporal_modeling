function E = entropy(logFPKM)

FPKM = 2.^logFPKM;
P = FPKM./sum(FPKM);
logP = log(P);
E = sum(P.*logP);


logFileName = 'sim2.txt';

logSize = 1000;
avgCallLength = 30;
callVar = 20;


a = round(rand(logSize,1)*logSize*100)
a = sort(a);
b = round((randn(logSize,1)*callVar)+avgCallLength);
b(b<5)=5;

calls = [a b];

fd = fopen(logFileName,'w');
fprintf(fd,'%d\t%d\n',calls');
fclose(fd);

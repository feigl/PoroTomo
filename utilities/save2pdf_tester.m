clear all
close all
histogram(randn([1000,1]));
%xlabel('X label can be very, very, very, very, very, very, very, very, very,very, very, very, very, very, very, very, very, very,very long');
xlabel('X label can be long');
ylabel('n');
save2pdf('test.pdf',gcf,1200);
%printpdf('test.pdf',1200);
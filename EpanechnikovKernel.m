function [ K ] = EpanechnikovKernel( dimension,inInds, x_metric )
%EpanechnikovKernel ���������ά����x�Ķ���ֵ����˺���ֵ
% K = ( (2.*pi).^(-dimension/2) ).*exp( -1/2.*(x_metric.^2) );
K = zeros(1,size(x_metric,2));
no_zero      = find(x_metric <= 1);              
K(no_zero) = 1./2 .*  (1./inInds) .* (dimension + 2) .* (1 - x_metric(no_zero));
% x_metric = x_metric ./ max(x_metric);
% no_zero      = find(x_metric <= 1);              
% K = 0.5 ./inInds .* (dimension + 2) .* (1 - x_metric);
end


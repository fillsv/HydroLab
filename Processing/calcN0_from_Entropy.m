function [N0] = calcN0_from_Entropy(signal)

% norm_signal = mat2gray(signal);
% En2 = entropy(norm_signal);
     signal =   signal(~isnan(signal));
     IIfII = sum(abs(signal).^2, 'all');
     %IIfII = sum(abs(omega), 'all')^2;
     p_i=(abs(signal).^2)/IIfII;
     ind = find(p_i ~= 0);
     
     En2=-1*sum(p_i(ind).*log(p_i(ind)), 'all');
     N0 = ceil(exp(En2));


     
%      keep4 = wentropy(C(:), 'shannon')/wentropy(C(:), 'norm', 2);
end
%@D
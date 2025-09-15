function action = chooseAction(Q, state, U_bins, h_bins, EbNo_bins)
    [~,iU] = min(abs(U_bins - state(1)));
    [~,iH] = min(abs(h_bins - state(2)));
    [~,iE] = min(abs(EbNo_bins - state(3)));
    [~,action] = max(Q(iU,iH,iE,:));
    action = action - 1; % 0 = SC, 1 = MRC
end

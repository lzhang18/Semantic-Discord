load demo_data

min_sigma = 0.005;

[D_BSF, i_BSF, j_BSF, l_BSF, ratio] = SemanticDiscord(ts, context_len, target_len,min_sigma);

[val, idx] = max(D_BSF);

disp(['distance value: ' num2str(val)])
disp(['context idx I: ' num2str(i_BSF(idx))])
disp(['context idx J: ' num2str(j_BSF(idx))])
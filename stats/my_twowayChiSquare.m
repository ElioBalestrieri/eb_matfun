function ChiStruct = my_twowayChiSquare(input_table)
%% two way table Chi Square Test

% convert table into array
input_arr = table2array(input_table);

% compute expected values
sum_rows = sum(input_arr); % along first dimension, sum all rows together.
sum_cols = sum(input_arr,2); % along second dimension, sum all cols together.
tot_occ = sum(sum(input_arr));
expected_vals = sum_cols * sum_rows / tot_occ;

ChiStat = sum(sum(((input_arr-expected_vals).^2)./expected_vals));

% get dof
chigof = prod(size(input_arr)-[1 1]);

P = 1-chi2cdf(ChiStat, chigof);

% append results (and input table to the final structure)
ChiStruct.TwoWaytable = input_table;
ChiStruct.p = P;
ChiStruct.ChiStat = ChiStat;
ChiStruct.ExpectedValues = expected_vals;


end


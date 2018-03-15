for a=1:numel(epochNames),
    relDays = cellfun(@(x) any(x==a),epochTypes);
    relDays = days(relDays);


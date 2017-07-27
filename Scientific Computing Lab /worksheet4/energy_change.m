function result = energy_change(solution1, solution2)

result = trace(solution2'*solution2) - trace(solution1'*solution1);
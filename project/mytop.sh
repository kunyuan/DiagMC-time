aa=$(ps axu | grep -E "run_loop|gamma3" | grep -v "sh" | grep -v "grep" | awk '{printf "%s,", $2}')
top -p ${aa:0:-1}

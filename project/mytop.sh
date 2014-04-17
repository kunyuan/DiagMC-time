aa=$(ps axu | grep -E "run_loop|gamma3" | grep -v "sh" | grep -v "grep" | awk '{printf "%s,", $2}')
top -pid ${aa:0:${#aa}-1}

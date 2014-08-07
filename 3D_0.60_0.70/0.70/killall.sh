aa=$(ps axu | grep -E "run_loop|gamma3" | grep -v "grep" | awk '{print $2}')
echo $aa
kill -9 $aa

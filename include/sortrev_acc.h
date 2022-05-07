  if(right == left)
    return;
  if(right - left > 16){
    l_hold = left + (int) ((right - left) * rand()/(1.*RAND_MAX));
    pivot = numbers[l_hold];
    ipivot = indices[l_hold];
    numbers[l_hold]=numbers[left];
    indices[l_hold]=indices[left];
    numbers[left] = pivot;
    indices[left] = ipivot;
    l_hold = left;
    r_hold = right;
    while (left < right){
      while ((numbers[right] <= pivot) && (left < right))
	right--;
      if (left != right){
	numbers[left] = numbers[right];
	indices[left] = indices[right];
	left++;}
      while ((numbers[left] >= pivot) && (left < right))
	left++;
      if (left != right){
	numbers[right] = numbers[left];
	indices[right] = indices[left];
	right--;}}
    numbers[left] = pivot;
    indices[left] = ipivot;
    right = (l_hold + r_hold)/2;
    if(left < right){
      while(numbers[left+1] == pivot && left < right)
	left++;}
    else{
      while(numbers[left-1] == pivot && left > right)
	left--;}      
  }
  else{
    pivot = numbers[left];
    ipivot = indices[left];
    l_hold = left;
    r_hold = right;
    while (left < right){
      while ((numbers[right] <= pivot) && (left < right))
	right--;
      if (left != right){
	numbers[left] = numbers[right];
	indices[left] = indices[right];
	left++;}
      while ((numbers[left] >= pivot) && (left < right))
	left++;
      if (left != right){
	numbers[right] = numbers[left];
	indices[right] = indices[left];
	right--;}}
    numbers[left] = pivot;
    indices[left] = ipivot;
  }

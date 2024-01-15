if(right == left)
  return;
l_hold = left;
r_hold = right;
if(right - left > 8){
  while(numbers[left]<numbers[left+1] && left< right)
    left++;
  if(left == right) return;
  while(numbers[right]>numbers[right-1] && right > left)
    right--;
  left = (left + right)/2;
  pivot = numbers[left];
  numbers[left]=numbers[l_hold];
  numbers[l_hold]=pivot;
  right = r_hold;
  left = l_hold;
  while (left < right){
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right){
      numbers[left] = numbers[right];
      left++;}
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right){
      numbers[right] = numbers[left];
      right--;}}
  numbers[left] = pivot;
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
   while (left < right){
     while ((numbers[right] >= pivot) && (left < right))
       right--;
     if (left != right){
       numbers[left] = numbers[right];
       left++;}
     while ((numbers[left] <= pivot) && (left < right))
       left++;
     if (left != right){
       numbers[right] = numbers[left];
       right--;}}
   numbers[left] = pivot;
 }

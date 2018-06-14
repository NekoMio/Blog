<?php
 
if(!isset($_GET['img']) || !file_exists(substr($_GET['img'],1))){
 header("HTTP/1.1 404 Not Found");
 echo "<h1>Not Found!</h1>";
 exit;
}
 
$img = substr($_GET['img'],1);
$img_webp = $img.'.webp';
if(!file_exists($img_webp)){
 `cwebp -q 100 $img -o $img_webp`;
}
 
header('Content-type: image/webp');
readfile($img_webp);
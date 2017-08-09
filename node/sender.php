<?php
if (isset($_GET['yl'])){
	$ts = time();
	$sig = genSig($ts);
	if (isset($_GET['er'])){
		$open = "http://localhost/Mining/Hpc/errnotify/id/{$_GET['yl']}/sig/{$sig}/t/{$ts}/err/{$_GET['er']}";
	}else{
		$open = "http://localhost/Mining/Hpc/done/id/{$_GET['yl']}/sig/{$sig}/t/{$ts}";
	}
	echo file_get_contents($open);
}else{
	echo '123456789';
}

function genSig($time) {
	$salt = '{SALT}';
	$toBuild = $time.$salt;
	 
	for ($i = 0; $i < 5; $i++) {
		$toBuild = md5($toBuild);
	}
	return $toBuild;
}
?>
<?php
if (isset($_GET['seq']) && isset($_GET['fb']) && isset($_GET['tamp']) && isset($_GET['sig'])){
	$ts = time();
	if (abs($ts - $_GET['tamp']) <= 600 && checkKey($ts, $_GET['sig'])){
		if ($_GET['seq'] == 1) {
			if (isset($_GET['url'])){
				passthru('python {RUNPATH}/dispatch.py 1 '.$_GET['url'].' '.$_GET['fb']);
			}else{
				echo '>Fail to check your format';
			}
		}elseif ($_GET['seq'] == 2){
			if (isset($_GET['ru'])){
				passthru('python {RUNPATH}/dispatch.py 2 '.$_GET['url'].' '.$_GET['fb']);
			}else{
				echo '>Fail to check your format';
			}
		}else{
			echo '>Fail to check your Ex type';
		}
	}else{
		echo '>Fail to check your signature.'
	}
}else{
	echo '>Are you attacking RNA Editing Plus server?';
}

function checkKey($timestamp, $key) {
	$salt = '{SALT}';
	$toBuild = $timestame.$salt;
	
	for ($i = 0; $i < 5; $i++) {
		$toBuild = md5($toBuild);
	}
	
	if ($key == $toBuild) {
		return 1;
	}else{
		return 0;
	}
}
?>
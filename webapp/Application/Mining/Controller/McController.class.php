<?php
namespace Mining\Controller;
use Think\Controller;
class McController extends Controller {
	/**
	 * Method 1
	 */
	public function vc($mail, $code) {
		$this->assign('mail', $mail);
		$this->assign('vcode', $code);
		$content = $this->fetch('Mc/vc');
		$m = think_send_mail($mail, $mail, 'Please active your query @ RNA Editing Plus', $content);
		if ($m === True){
			return 1;
		}else{
			return $m;
		}
	}
	
	/**
	 * Method 2
	 */
	public function d($mail, $code) {
		$this->assign('mail', $mail);
		$this->assign('vcode', $code);
		$content = $this->fetch('Mc/d');
		$m = think_send_mail($mail, $mail, 'Congratulations! Your job is finished (RNA Editing Plus)', $content);
		if ($m === True){
			return 1;
		}else{
			return $m;
		}
	}

	public function downerr($mail, $code) {
		$this->assign('mail', $mail);
		$this->assign('vcode', $code);
		$content = $this->fetch('Mc/derr');
		$m = think_send_mail($mail, $mail, 'Failed to check your inputs (RNA Editing Plus)', $content);
		if ($m === True){
			return 1;
		}else{
			return $m;
		}
	}

	public function operr($mail, $code) {
		$this->assign('mail', $mail);
		$this->assign('vcode', $code);
		$content = $this->fetch('Mc/err');
		$m = think_send_mail($mail, $mail, 'Failed to analyse your samples (RNA Editing Plus)', $content);
		if ($m === True){
			return 1;
		}else{
			return $m;
		}
	}

	public function queueMail($type, $mail, $code){
		$mailQ = M('mail_queue');
		$data['type'] = $type;
		$data['mail'] = $mail;
		$data['vcode'] = $code;
		$data['status'] = 0;
		if ($mailQ->data($data)->add()){
			return 1;
		}else{
			$this->error('Error when pushing task into our queue.');
		}
	}

	public function sendMail(){
		
		set_time_limit(600);
		if($_SERVER['HTTP_HOST'] == 'localhost'){
			$driver = M('mail_queue');
			$list = $driver->where('status=0 AND times < 3')->limit(10)->select();
			if (count($list) > 0){
				foreach ($list as $sub) {
					if ($sub['type'] == 1){
						$ret = $this->vc($sub['mail'], $sub['vcode']);
					}elseif ($sub['type'] == 2){
						$ret = $this->d($sub['mail'], $sub['vcode']);
					}elseif ($sub['type'] == 3){
						$ret = $this->downerr($sub['mail'], $sub['vcode']);
					}elseif ($sub['type'] == 4){
						$ret = $this->operr($sub['mail'], $sub['vcode']);
					}elseif ($sub['type'] == -1){
						$ret = $this->operr($sub['mail'], $sub['vcode']);
					}
					if ($ret != 1){
						$data['times'] = $sub['times'] + 1;
						$data['msg'] = $ret;
					}else{
						$data['status'] = 1;
					}
					$driver->where('id='.$sub['id'])->save($data);
				}
			}
		}
	}
}
?>
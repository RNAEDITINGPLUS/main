<?php
namespace Mining\Controller;
use Think\Controller;
class HpcController extends Controller {
    public function transfer() {
    	$job = M('jobs');
    	$ts = time();
    	if ($job->autoCheckToken($_POST) && check_verify(I('verify'))) {
    		$fbTag = I('fb') == 1 ? '&fb=1' : '&fb=0';
    		$data['user'] = I('user');
    		$data['st'] = I('seq');
    		$data['fb'] = I('fb');
    		$data['tissue'] = I('tissue');
			$data['vcode'] = create_guid();
    		$data['description'] = I('description');
			if ($_FILES['es']['name']){
				$upload = new \Think\Upload();
				$upload->maxSize   =     10485760;
				$upload->exts      =     array('txt', 'gtf', 'csv');
				$upload->rootPath  =      './Inner/';
				$upload->savePath  =      '';
				$info   =   $upload->uploadOne($_FILES['es']);
				if(!$info) {
					$this->error($upload->getError());
				}else{
					if ($_FILES['expr']['name']){
						$eInfo = $upload->uploadOne($_FILES['exprexpr']);
						if(!$eInfo){
							$this->error($upload->getError());
						}else{
							$data['es'] = 'http://host.rnaeditplus.org/Inner/'.$info['savepath'].$info['savename'];
							$data['expr'] = 'http://host.rnaeditplus.org/Inner/'.$eInfo['savepath'].$eInfo['savename'];
						}
					}else{
						$data['es'] = 'http://host.rnaeditplus.org/Inner/'.$info['savepath'].$info['savename'];
					}
					
				}
				//$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?es='.$data['es'].$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
			}else{
				if (I('du') && I('ru')){
    				$data['job'] = I('du').';'.I('ru');
    				//$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?seq='.I('seq').'&a='.I('du').'&b='.I('ru').$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
    			}else if(I('ru')){
    				$data['job'] = I('ru');
    				//$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?seq='.I('seq').'&b='.I('ru').$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
    			}
			}
    		//$ct = rtrim(file_get_contents($innerListener));
    		//if (is_numeric($ct)) {
    		//	$data['trace'] = $ct;
    		$job->add($data);
			$notify = A('Mc');
			$n = $notify->queueMail(1, $data['user'], $data['vcode']);
			if($n == 1){
				$this->assign('mail', $data['user']);
				$this->display();
			}else{
				$this->error($n);
			}
    		//	$this->success('Your task had beed pushed into our queue.');
    		//}else{
    		//	$this->error($ct);
    		//}
    	}else{
    		$this->error("> Are you attacking RNA Editing Plus server?<br />This form is undes Bro Li's protection.");
    	}
    }
    public function proceed($vc){
		$job = M('jobs');
		$map = array(
			'vcode'	=>	$vc,
			'status'	=>	0,
		);
		$ret = $job->where($map)->find();
		if ($ret){
			
			$fbTag = $ret['fb'] == 1 ? '&fb=1' : '&fb=0';
			$ts = time();
			
			if ($ret['es'] != ''){
				if ($ret['expr'] != ''){
					$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?es='.$ret['es'].'&ex='.$ret['expr'].$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
				}else{
					$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?es='.$ret['es'].$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
				}
				
			}elseif(strpos($ret['job'], ';')){
				$u = explode(';', $ret['job']);
				$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?seq=2&a='.$u[0].'&b='.$u[1].$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
			}elseif($ret['job']){
				$innerListener = 'http://h.rnaeditplus.org/dispatch/autorun/listener.php?seq='.I('seq').'&b='.$ret['job'].$fbTag.'&tamp='.$ts.'&sig='.$this->genSig($ts);
			}
			$ct = rtrim(file_get_contents($innerListener));
    		if (is_numeric($ct)) {
				$map['id'] = $ret['id'];
    			$data['trace'] = $ct;
				$data['status'] = 1;
				$job->where($map)->save($data);
				$this->success('Your task had beed pushed into our queue.', '/Viewer/Index/index');
			}else{
				$this->error($ct);
			}
			
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is undes Bro Li's protection.");
		}
	}
	
	public function done($id, $sig, $t){
		if ($this->checkKey($t, $sig)){
			$job = M('jobs');
			$map = array('trace'	=>	$id);
			$user = $job->where($map)->getField('id, user, vcode');
			$notify = A('Mc');
			$user = array_pop($user);
			$n = $notify->queueMail(2, $user["user"], $user["vcode"]);
			if($n == 1){
				echo 0;
			}else{
				echo $n;
			}
		}else{
			echo "Sig error;";
			#$this->error("> Are you attacking RNA Editing Plus server?<br />This form is undes Bro Li's protection.");
		}
	}
	public function errnotify($id, $sig, $t, $err){
		if ($this->checkKey($t, $sig)){
			$job = M('jobs');
			$map = array('trace'	=>	$id);
			$user = $job->where($map)->getField('id, user, vcode');
			$notify = A('Mc');
			$user = array_pop($user);
			$n = $notify->queueMail($err, $user["user"], $user["vcode"]);
			if($n == 1){
				echo 0;
			}else{
				echo $n;
			}
		}else{
			echo "Sig error;";
			#$this->error("> Are you attacking RNA Editing Plus server?<br />This form is undes Bro Li's protection.");
		}
	}
    private function genSig($time) {
    	$salt = 'Q#77q0$&xyp18o8%';
    	$toBuild = $time.$salt;
    	
    	for ($i = 0; $i < 5; $i++) {
    		$toBuild = md5($toBuild);
    	}
    	return $toBuild;
    }
    private function checkKey($timestamp, $key) {
    	$salt = 'Q#77q0$&xyp18o8%';
    	$toBuild = $timestamp.$salt;
    
    	for ($i = 0; $i < 5; $i++) {
    		$toBuild = md5($toBuild);
    	}
    	if ($key == $toBuild) {
    		return 1;
    	}else{
    		return 0;
    	}
    }
	public function testRed(){
		$this->success('Your task had beed pushed into our queue.', 'Viewer/Index/index');
	}
}
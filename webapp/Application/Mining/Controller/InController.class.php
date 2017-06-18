<?php
namespace Mining\Controller;
use Think\Controller;
class InController extends Controller {
    public function form(){
		$this->assign('tisl', $this->getTissue());
        $this->display();
    }
    
    public function transfer() {
    	$job = M('innovation');
    	$ts = time();
    	if ($job->autoCheckToken($_POST) && check_verify(I('verify'))) {
    		$data['user'] = session('user');
    		$data['description'] = I('description');
    		if (I('du') && I('ru')){
    			$data['job'] = I('du').';'.I('ru');
    			$innerListener = 'http://10.117.228.126/dispatch/autorun/y/listener.php?a='.I('du').'&b='.I('ru').'&tamp='.$ts.'&sig='.$this->genSig($ts);
    		}else if(I('ru')){
    			$data['job'] = I('ru');
    			$innerListener = 'http://10.117.228.126/dispatch/autorun/y/listener.php?a='.I('ru').'&tamp='.$ts.'&sig='.$this->genSig($ts);
    		}
    		$ct = rtrim(file_get_contents($innerListener));
    		if (is_numeric($ct)) {
    			$data['jobid'] = $ct;
				$data['tissue'] = I('tissue');
    			$job->add($data);
    			$this->success('Your task had beed pushed into our queue.');
    		}else{
    			$this->error($ct);
    		}
    	}else{
    		print "> Are you attacking RNA Editing Plus server?<br />";
    		print "This form is undes Bro Li's protection.";
    	}
    }
    public function nhpc(){
		$job = M('innovation_nefuc');
		$data = I();
		echo $job->add($data);
	}
    private function genSig($time) {
    	$salt = 'Q#77q0$&xyp18o8%';
    	$toBuild = $time.$salt;
    	
    	for ($i = 0; $i < 5; $i++) {
    		$toBuild = md5($toBuild);
    	}
    	return $toBuild;
    }
    
    public function verify(){
    	$verify = new \Think\Verify();
    	$verify->useImgBg = true;
    	$verify->entry(1);
    }
	
	public function getJob(){
		$job = M('innovation_nefuc');
		$check=count($job->where('status!=0 AND status!=6')->select());
		if ($check < 3){
			$newJob = $job->where("status='0'")->order('id')->limit(1)->select();
			if (count($newJob) == 1){
				$newJobs = $job->where("id={$newJob[0]['id']}")->setField('status', 1);
				echo $newJob[0]['id']."\t".$newJob[0]['job'];
			}else{
				echo 'nonono';
			}
		}
		
	}
	
	public function updateJob($job, $field, $st){
		$j = M('innovation_nefuc');
		$newJobs = $j->where("id={$job}")->setField($field, $st);
	}
	
	private function getTissue(){
		$ti = M('tissue');
		return $ti->cache('tissue')->where('status=1')->select();
	}
}
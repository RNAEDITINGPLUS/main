<?php
namespace Mining\Controller;
use Think\Controller;
class IndexController extends Controller {
    public function index(){
		$this->assign('tisl', $this->getTissue());
        $this->display();
    }
	
	public function verify(){
    	$verify = new \Think\Verify();
    	$verify->useImgBg = true;
    	$verify->entry(1);
    }
	
	private function getTissue(){
		$ti = M('tissue');
		return $ti->cache('tissue')->where('status=1')->select();
	}
}
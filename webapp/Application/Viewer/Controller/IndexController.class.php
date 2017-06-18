<?php
namespace Viewer\Controller;
use Think\Controller;
class IndexController extends Controller {
    public function index(){
        $this->display();
    }
	
	public function verify(){
    	$verify = new \Think\Verify();
    	$verify->useImgBg = true;
    	$verify->entry(1);
    }
    public function datasets(){
    	$this->display();
    }
}
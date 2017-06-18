<?php
// +----------------------------------------------------------------------
// | ThinkPHP [ WE CAN DO IT JUST THINK IT ]
// +----------------------------------------------------------------------
// | Copyright (c) 2006-2014 http://thinkphp.cn All rights reserved.
// +----------------------------------------------------------------------
// | Licensed ( http://www.apache.org/licenses/LICENSE-2.0 )
// +----------------------------------------------------------------------
// | Author: liu21st <liu21st@gmail.com>
// +----------------------------------------------------------------------

defined('THINK_PATH') or exit();
return  array(
    /* 应用设定 */
    'APP_USE_NAMESPACE'     =>  true,
    'APP_SUB_DOMAIN_DEPLOY' =>  false,
    'APP_SUB_DOMAIN_RULES'  =>  array(),
    'APP_DOMAIN_SUFFIX'     =>  '',
    'ACTION_SUFFIX'         =>  '',
    'MULTI_MODULE'          =>  true,
    'MODULE_DENY_LIST'      =>  array('Common','Runtime'),
    'CONTROLLER_LEVEL'      =>  1,
    'APP_AUTOLOAD_LAYER'    =>  'Controller,Model',
    'APP_AUTOLOAD_PATH'     =>  '',

    /* Cookie设置 */
    'COOKIE_EXPIRE'         =>  0,
    'COOKIE_DOMAIN'         =>  '',
    'COOKIE_PATH'           =>  '/',
    'COOKIE_PREFIX'         =>  '',
    'COOKIE_SECURE'         =>  false,
    'COOKIE_HTTPONLY'       =>  '',

    /* 默认设定 */
    'DEFAULT_M_LAYER'       =>  'Model',
    'DEFAULT_C_LAYER'       =>  'Controller',
    'DEFAULT_V_LAYER'       =>  'View',
    'DEFAULT_LANG'          =>  'zh-cn',
    'DEFAULT_THEME'         =>  '',
    'DEFAULT_MODULE'        =>  'Viewer',
    'DEFAULT_CONTROLLER'    =>  'Index',
    'DEFAULT_ACTION'        =>  'index',
    'DEFAULT_CHARSET'       =>  'utf-8',
    'DEFAULT_TIMEZONE'      =>  'PRC',
    'DEFAULT_AJAX_RETURN'   =>  'JSON',
    'DEFAULT_JSONP_HANDLER' =>  'jsonpReturn',
    'DEFAULT_FILTER'        =>  'htmlspecialchars',

    /* 数据库设置 */
    'DB_TYPE'               =>  'mysql',     // Database type
    'DB_HOST'               =>  '', // Database host
    'DB_NAME'               =>  '',          // Database name
    'DB_USER'               =>  '',      // User name
    'DB_PWD'                =>  '',          // Password
    'DB_PORT'               =>  '',        // Port
    'DB_PREFIX'             =>  '',    // Prefix
    'DB_PARAMS'          	=>  array(),    
    'DB_DEBUG'  			=>  TRUE,
    'DB_FIELDS_CACHE'       =>  true,
    'DB_CHARSET'            =>  'utf8',
    'DB_DEPLOY_TYPE'        =>  0,
    'DB_RW_SEPARATE'        =>  false,
    'DB_MASTER_NUM'         =>  1,
    'DB_SLAVE_NO'           =>  '',

    /* 数据缓存设置 */
    'DATA_CACHE_TIME'       =>  0,
    'DATA_CACHE_COMPRESS'   =>  false,
    'DATA_CACHE_CHECK'      =>  false,
    'DATA_CACHE_PREFIX'     =>  '',
    'DATA_CACHE_TYPE'       =>  'File',
    'DATA_CACHE_PATH'       =>  TEMP_PATH,
    'DATA_CACHE_KEY'        =>  '',
    'DATA_CACHE_SUBDIR'     =>  false,
    'DATA_PATH_LEVEL'       =>  1,

    /* 错误设置 */
    'ERROR_MESSAGE'         =>  'Something wrong...',
    'ERROR_PAGE'            =>  '',
    'SHOW_ERROR_MSG'        =>  false,
    'TRACE_MAX_RECORD'      =>  100,

    /* 日志设置 */
    'LOG_RECORD'            =>  false,
    'LOG_TYPE'              =>  'File',
    'LOG_LEVEL'             =>  'EMERG,ALERT,CRIT,ERR',
    'LOG_FILE_SIZE'         =>  2097152,
    'LOG_EXCEPTION_RECORD'  =>  false,

    /* SESSION设置 */
    'SESSION_AUTO_START'    =>  true,    // 是否自动开启Session
    'SESSION_OPTIONS'       =>  array(), // session 配置数组 支持type name id path expire domain 等参数
    'SESSION_TYPE'          =>  '', // session hander类型 默认无需设置 除非扩展了session hander驱动
    'SESSION_PREFIX'        =>  '', // session 前缀
    //'VAR_SESSION_ID'      =>  'session_id',     //sessionID的提交变量

    /* 模板引擎设置 */
    'TMPL_CONTENT_TYPE'     =>  'text/html', // 默认模板输出类型
    'TMPL_ACTION_ERROR'     =>  THINK_PATH.'Tpl/dispatch_jump.tpl', // 默认错误跳转对应的模板文件
    'TMPL_ACTION_SUCCESS'   =>  THINK_PATH.'Tpl/dispatch_jump.tpl', // 默认成功跳转对应的模板文件
    'TMPL_EXCEPTION_FILE'   =>  THINK_PATH.'Tpl/think_exception.tpl',// 异常页面的模板文件
    'TMPL_DETECT_THEME'     =>  false,       // 自动侦测模板主题
    'TMPL_TEMPLATE_SUFFIX'  =>  '.html',     // 默认模板文件后缀
    'TMPL_FILE_DEPR'        =>  '/', //模板文件CONTROLLER_NAME与ACTION_NAME之间的分割符
    // 布局设置
    'TMPL_ENGINE_TYPE'      =>  'Think',     // 默认模板引擎 以下设置仅对使用Think模板引擎有效
    'TMPL_CACHFILE_SUFFIX'  =>  '.php',      // 默认模板缓存后缀
    'TMPL_DENY_FUNC_LIST'   =>  'echo,exit',    // 模板引擎禁用函数
    'TMPL_DENY_PHP'         =>  false, // 默认模板引擎是否禁用PHP原生代码
    'TMPL_L_DELIM'          =>  '{',            // 模板引擎普通标签开始标记
    'TMPL_R_DELIM'          =>  '}',            // 模板引擎普通标签结束标记
    'TMPL_VAR_IDENTIFY'     =>  'array',     // 模板变量识别。留空自动判断,参数为'obj'则表示对象
    'TMPL_STRIP_SPACE'      =>  true,       // 是否去除模板文件里面的html空格与换行
    'TMPL_CACHE_ON'         =>  true,        // 是否开启模板编译缓存,设为false则每次都会重新编译
    'TMPL_CACHE_PREFIX'     =>  '',         // 模板缓存前缀标识，可以动态改变
    'TMPL_CACHE_TIME'       =>  0,         // 模板缓存有效期 0 为永久，(以数字为值，单位:秒)
    'TMPL_LAYOUT_ITEM'      =>  '{__CONTENT__}', // 布局模板的内容替换标识
    'LAYOUT_ON'             =>  false, // 是否启用布局
    'LAYOUT_NAME'           =>  'layout', // 当前布局名称 默认为layout

    // Think模板引擎标签库相关设定
    'TAGLIB_BEGIN'          =>  '<',  // 标签库标签开始标记
    'TAGLIB_END'            =>  '>',  // 标签库标签结束标记
    'TAGLIB_LOAD'           =>  true, // 是否使用内置标签库之外的其它标签库，默认自动检测
    'TAGLIB_BUILD_IN'       =>  'cx', // 内置标签库名称(标签使用不必指定标签库名称),以逗号分隔 注意解析顺序
    'TAGLIB_PRE_LOAD'       =>  '',   // 需要额外加载的标签库(须指定标签库名称)，多个以逗号分隔 
    
    /* URL设置 */
    'URL_CASE_INSENSITIVE'  =>  true,   // 默认false 表示URL区分大小写 true则表示不区分大小写
    'URL_MODEL'             =>  1,       // URL访问模式,可选参数0、1、2、3,代表以下四种模式：
    // 0 (普通模式); 1 (PATHINFO 模式); 2 (REWRITE  模式); 3 (兼容模式)  默认为PATHINFO 模式
    'URL_PATHINFO_DEPR'     =>  '/',	// PATHINFO模式下，各参数之间的分割符号
    'URL_PATHINFO_FETCH'    =>  'ORIG_PATH_INFO,REDIRECT_PATH_INFO,REDIRECT_URL', // 用于兼容判断PATH_INFO 参数的SERVER替代变量列表
    'URL_REQUEST_URI'       =>  'REQUEST_URI', // 获取当前页面地址的系统变量 默认为REQUEST_URI
    'URL_HTML_SUFFIX'       =>  'html',  // URL伪静态后缀设置
    'URL_DENY_SUFFIX'       =>  'ico|png|gif|jpg', // URL禁止访问的后缀设置
    'URL_PARAMS_BIND'       =>  true, // URL变量绑定到Action方法参数
    'URL_PARAMS_BIND_TYPE'  =>  0, // URL变量绑定的类型 0 按变量名绑定 1 按变量顺序绑定
    'URL_PARAMS_FILTER'     =>  false, // URL变量绑定过滤
    'URL_PARAMS_FILTER_TYPE'=>  '', // URL变量绑定过滤方法 如果为空 调用DEFAULT_FILTER
    'URL_ROUTER_ON'         =>  false,   // 是否开启URL路由
    'URL_ROUTE_RULES'       =>  array(), // 默认路由规则 针对模块
    'URL_MAP_RULES'         =>  array(), // URL映射定义规则

    /* 系统变量名称设置 */
    'VAR_MODULE'            =>  'm',     // 默认模块获取变量
    'VAR_ADDON'             =>  'addon',     // 默认的插件控制器命名空间变量
    'VAR_CONTROLLER'        =>  'c',    // 默认控制器获取变量
    'VAR_ACTION'            =>  'a',    // 默认操作获取变量
    'VAR_AJAX_SUBMIT'       =>  'ajax',  // 默认的AJAX提交变量
    'VAR_JSONP_HANDLER'     =>  'callback',
    'VAR_PATHINFO'          =>  's',    // 兼容模式PATHINFO获取变量例如 ?s=/module/action/id/1 后面的参数取决于URL_PATHINFO_DEPR
    'VAR_TEMPLATE'          =>  't',    // 默认模板切换变量
    'VAR_AUTO_STRING'		=>	false,	// 输入变量是否自动强制转换为字符串 如果开启则数组变量需要手动传入变量修饰符获取变量

    'HTTP_CACHE_CONTROL'    =>  'private',  // 网页缓存控制
    'CHECK_APP_DIR'         =>  true,       // 是否检查应用目录是否创建
    'FILE_UPLOAD_TYPE'      =>  'Local',    // 文件上传方式
    'DATA_CRYPT_TYPE'       =>  'Think',    // 数据加密方式

);

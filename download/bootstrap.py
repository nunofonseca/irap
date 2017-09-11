<!doctype html>
<!-- Server: sfs-consume-17 -->

<!--[if lt IE 7 ]> <html lang="en" class="no-js ie6"> <![endif]-->
<!--[if IE 7 ]>    <html lang="en" class="no-js ie7"> <![endif]-->
<!--[if IE 8 ]>    <html lang="en" class="no-js ie8"> <![endif]-->
<!--[if IE 9 ]>    <html lang="en" class="no-js ie9"> <![endif]-->
<!--[if (gt IE 9)|!(IE)]>--> <html lang="en" class="no-js"> <!--<![endif]-->
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
        
        
        <script type='text/javascript'>
            /*global unescape, window, SF*/
            // Setup our namespace
            if (!window.SF) { window.SF = {}; }
            if (!window.net) { window.net = {}; }
            if (!window.net.sf) { window.net.sf = {}; }
            SF.Ads = {};
            SF.cdn = '//a.fsdn.com/con';
            SF.deploy_time = '1504810730';
            
            SF.Breakpoints = {
              small: 0,
              medium: 640,
              leaderboard: 743,
              billboard: 985,
              large: 1053,
              xlarge: 1295,
              xxlarge: 1366
            };
            SF.initial_breakpoints_visible = {};
            for (var bp in SF.Breakpoints) {
                if (!SF.Breakpoints.hasOwnProperty(bp)) {
                    continue;
                }
                SF.initial_breakpoints_visible[bp] = !window.matchMedia || window.matchMedia('(min-width: ' + SF.Breakpoints[bp] + 'px)').matches;
            }
        </script><script type='text/javascript'>
            SF.Ads.prebidOptions = {
                showIndicators: false,
                analytics: true,
                timeout: 650,
                timeouts_by_bids: {650: 3, 2000: 1, 3000: 0},
                };

            SF.Ads.prebidUnits = [];
            if (SF.initial_breakpoints_visible.leaderboard) {
                SF.Ads.prebidUnits.push({
                    bids: [{'params': {'tagid': '364648'}, 'bidder': 'sovrn'},
                            {'params': {'cp': 558092, 'cf': '728x90', 'ct': '472223'}, 'bidder': 'pulsepoint'},
                            {'params': {'placementId': 8841121}, 'bidder': 'brealtime'},
                            {'params': {'siteId': '103240', 'position': 'atf', 'sizes': [2], 'accountId': '15680', 'zoneId': '486110'}, 'bidder': 'rubicon'},
                            {'params': {'placementId': 9265078}, 'bidder': 'appnexus'},
                            {'params': {'placement': 4224498, 'network': '10676.1'}, 'bidder': 'aol'},
                            {'params': {'siteID': '188656', 'id': '7'}, 'bidder': 'indexExchange'},
                            {'params': {'placementId': '65725'}, 'bidder': 'rhythmone'},
                            ],
                    code: 'div-gpt-ad-1393435113147-0',
                    tag: 'SF_ProjectFiles_728x90_A',
                    
                    sizes: [728, 90]
                });
            }
            if (true) {
                SF.Ads.prebidUnits.push({
                    bids: [{'params': {'tagid': '364646'}, 'bidder': 'sovrn'},
                            {'params': {'cp': 558092, 'cf': '300x250', 'ct': '472221'}, 'bidder': 'pulsepoint'},
                            {'params': {'placementId': 9319285}, 'bidder': 'brealtime'},
                            {'params': {'siteId': '103240', 'position': 'atf', 'sizes': [15, 10], 'accountId': '15680', 'zoneId': '486110'}, 'bidder': 'rubicon'},
                            {'params': {'placementId': 9265076}, 'bidder': 'appnexus'},
                            {'params': {'placement': 4224505, 'network': '10676.1'}, 'bidder': 'aol'},
                            {'params': {'siteID': '188654', 'id': '5'}, 'bidder': 'indexExchange'},
                            {'params': {'placementId': '65725'}, 'bidder': 'rhythmone'},
                            ],
                    code: 'div-gpt-ad-1392147725721-0',
                    tag: 'SF_ProjectFiles_300x250_A',
                    
                    sizes: [[300, 250], [300, 600], [300, 1050]]
                });
            }
            if (true) {
                SF.Ads.prebidUnits.push({
                    bids: [{'params': {'tagid': '364647'}, 'bidder': 'sovrn'},
                            {'params': {'cp': 558092, 'cf': '300x250', 'ct': '472222'}, 'bidder': 'pulsepoint'},
                            {'params': {'placementId': 8829468}, 'bidder': 'brealtime'},
                            {'params': {'siteId': '103240', 'sizes': [15], 'accountId': '15680', 'zoneId': '486112'}, 'bidder': 'rubicon'},
                            {'params': {'placementId': 9265082}, 'bidder': 'appnexus'},
                            {'params': {'placement': 4224503, 'network': '10676.1'}, 'bidder': 'aol'},
                            {'params': {'siteID': '188655', 'id': '6'}, 'bidder': 'indexExchange'},
                            {'params': {'placementId': '65725'}, 'bidder': 'rhythmone'},
                            ],
                    code: 'div-gpt-ad-1392148208789-0',
                    tag: 'SF_ProjectFiles_300x250_B',
                    
                    sizes: [300, 250]
                });
            }
            SF.Ads.prebidAdjustments = {"bidder_deflations": {"komoona": 0.92, "indexexchange": 0.9, "sovrn": 0.9, "aardvark": 0.99, "aol": 0.95, "brealtime": 0.88, "rubiconlite": 0.9, "pulsepoint": 0.89, "rubicon": 0.98, "indexex#hange": 0.99, "springserve": 0.001, "appnexus": 0.97, "rhythmone": 0.95}, "inflation": 1.2, "floor": 0.02};
        </script>
        <script type="text/javascript" id="pbjs_script" data-dom="https://d3tglifpd8whs6.cloudfront.net"  src="//a.fsdn.com/con/js/sftheme/vendor/bizx-prebid.js?1504810730" ></script>
        
        <meta id="project_name" name="project_name" content="fusioncatcher">
        
        <meta name="description" content="Somatic fusion-genes finder for RNA-seq data">
        <meta name="keywords" content="Bio-Informatics,  Open Source, Open Source Software, Development, Community, Source Code, Secure,  Downloads, Free Software">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>FusionCatcher -  Browse Files at SourceForge.net</title>
        <link rel="shortcut icon" href="//a.fsdn.com/con/img/sftheme/favicon.ico">
        
        <script type="text/javascript">
            /*global unescape, window, console, jQuery, $, net, SF, DD_belatedPNG, ga */
            if (!window.SF) {
                window.SF = {};
            }
        </script>

        <script type="text/javascript">
            SF.EU_country_codes = ["BE", "FR", "BG", "DK", "VG", "WF", "HR", "BM", "DE", "HU", "JE", "FI", "FK", "YT", "NL", "PT", "CW", "NC", "LV", "RE", "LT", "LU", "PF", "GI", "TF", "RO", "PN", "TC", "PL", "PM", "GS", "GR", "GP", "EE", "IT", "GG", "CZ", "CY", "SX", "IO", "AT", "AW", "AX", "GL", "IE", "KY", "ES", "ME", "MF", "BL", "GF", "SK", "MT", "SI", "SH", "MQ", "MS", "AI", "SE", "GB"];
            SF.unknown_country_codes = ["", "A1", "A2", "O1"];
        </script>
        
<script src="//a.fsdn.com/con/js/sftheme/vendor/modernizr.3.3.1.custom.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/vendor/jquery-1.11.1.min.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.core.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/bootstrap.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/sticky.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/shared_head.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme-typescript/compliance.js?1504810730" type="text/javascript"></script>


<!--[if lt IE 7 ]>
  <script src="//a.fsdn.com/con/js/sftheme/vendor/dd_belatedpng.js"></script>
  <script> DD_belatedPNG.fix('img, .png_bg'); //fix any <img> or .png_bg background-images </script>
<![endif]-->


<link href='//fonts.googleapis.com/css?family=Ubuntu:regular' rel='stylesheet' type='text/css'>
<style type="text/css">
    @font-face {
        font-family: "Pictos";
        src: url('//a.fsdn.com/con/css/fonts/sftheme/pictos-web.eot');
        src: local("☺"), url('//a.fsdn.com/con/css/fonts/sftheme/pictos-web.woff') format('woff'), url('//a.fsdn.com/con/css/fonts/sftheme/pictos-web.ttf') format('truetype'), url('//a.fsdn.com/con/css/fonts/sftheme/pictos-web.svg') format('svg');
    }
</style>
 
<link rel="stylesheet" href="//a.fsdn.com/con/css/sf-files.css?1504810730" type="text/css">



        <style type="text/css">.H9effeddc4f6fef7a225a171e663950d2074984de { display: none !important; }</style>

        
<script type="text/javascript">
    SF.adblock = true;
</script><script type='text/javascript'>
       /*global Foundation */
           /*global googletag, bizxPrebid */
       var gptadslots=[];
       var gptadHandlers={};
       var gptadRenderers=[];
       var gptadComplements={};

       googletag.cmd.push(function() {
           var leaderboard = googletag.sizeMapping()
               .addSize([970, 200], [[970, 250], [728, 90]])
               .addSize([728, 200], [[728, 90]])
               .build();
           var leaderboardInContent = googletag.sizeMapping()
               .addSize([1280, 200], [[970, 250], [728, 90]])
               .addSize([728, 200], [[728, 90]])
               .build();

            //prebid_log('GPT push define slots and targeting');
            gptadslots.push(googletag.defineSlot('/41014381/Sourceforge/SF_ProjectFiles_728x90_A',[728, 90],'div-gpt-ad-1393435113147-0')
                                                    .addService(googletag.pubads())
                                                        .setTargeting('tpc', ['fusioncatcher', 'python', 'scientific', 'bioinformatics'])
                                                        .setTargeting('shortname', 'fusioncatcher')
                                                        .setTargeting('aud', ['scienceresearch'])
                                                        .setTargeting('mirror', 'False')
                                                        .setTargeting('dc_ref', 'https://sourceforge.net/projects/fusioncatcher/files/')
                                                        .setTargeting('sz', '728x90')
                                                        .setTargeting('page_type', 'pg_files'));
            gptadslots.push(googletag.defineSlot('/41014381/Sourceforge/SF_ProjectFiles_300x250_A',[[300, 250], [300, 600], [300, 1050]],'div-gpt-ad-1392147725721-0')
                                                    .addService(googletag.pubads())
                                                        .setTargeting('tpc', ['fusioncatcher', 'python', 'scientific', 'bioinformatics'])
                                                        .setTargeting('shortname', 'fusioncatcher')
                                                        .setTargeting('aud', ['scienceresearch'])
                                                        .setTargeting('mirror', 'False')
                                                        .setTargeting('dc_ref', 'https://sourceforge.net/projects/fusioncatcher/files/')
                                                        .setTargeting('sz', '300x250,300x600,300x1050')
                                                        .setTargeting('page_type', 'pg_files'));
            gptadslots.push(googletag.defineSlot('/41014381/Sourceforge/SF_ProjectFiles_300x250_B',[300, 250],'div-gpt-ad-1392148208789-0')
                                                    .addService(googletag.pubads())
                                                        .setTargeting('tpc', ['fusioncatcher', 'python', 'scientific', 'bioinformatics'])
                                                        .setTargeting('shortname', 'fusioncatcher')
                                                        .setTargeting('aud', ['scienceresearch'])
                                                        .setTargeting('mirror', 'False')
                                                        .setTargeting('dc_ref', 'https://sourceforge.net/projects/fusioncatcher/files/')
                                                        .setTargeting('sz', '300x250')
                                                        .setTargeting('page_type', 'pg_files'));
            gptadslots.push(googletag.defineSlot('/41014381/Sourceforge/SF_ProjectFiles_HubIcon_200x90_A',[200, 90],'div-gpt-ad-1392148098424-0')
                                                    .addService(googletag.pubads())
                                                        .setTargeting('tpc', ['fusioncatcher', 'python', 'scientific', 'bioinformatics'])
                                                        .setTargeting('shortname', 'fusioncatcher')
                                                        .setTargeting('aud', ['scienceresearch'])
                                                        .setTargeting('mirror', 'False')
                                                        .setTargeting('dc_ref', 'https://sourceforge.net/projects/fusioncatcher/files/')
                                                        .setTargeting('sz', '200x90')
                                                        .setTargeting('page_type', 'pg_files'));
            gptadslots.push(googletag.defineSlot('/7346874/SF-300x250',[300, 250],'div-gpt-ad-1392148208795-0')
                                                    .addService(googletag.pubads())
                                                        .setTargeting('tpc', ['fusioncatcher', 'python', 'scientific', 'bioinformatics'])
                                                        .setTargeting('shortname', 'fusioncatcher')
                                                        .setTargeting('aud', ['scienceresearch'])
                                                        .setTargeting('mirror', 'False')
                                                        .setTargeting('dc_ref', 'https://sourceforge.net/projects/fusioncatcher/files/')
                                                        .setTargeting('sz', '300x250')
                                                        .setTargeting('page_type', 'pg_files'));
            

            googletag.pubads().setTargeting('requestSource', 'GPT');
            googletag.pubads().enableAsyncRendering();

            googletag.pubads().collapseEmptyDivs();
            googletag.pubads().addEventListener('slotRenderEnded', function(event) {
                var unitName = event.slot.getName();
                if ( unitName in gptadHandlers ) {
                   for (var i = 0; i < gptadHandlers[unitName].length; i++) {
                       try {
                           SF.Ads.RenderHandlers[gptadHandlers[unitName][i]].call(this, event);
                       } catch (e) {
                       }
                   }
                }
            });
            googletag.pubads().addEventListener('impressionViewable', SF.Ads.RenderHandlers.viewabilityInstrumentation);
            
            googletag.pubads().addEventListener('slotRenderEnded', SF.Ads.listenerForBlockThis);
            bizxPrebid.Ads.pushToGoogle();

            googletag.enableServices();
        });
   </script> 

        

        
<!-- CCM Tag -->
<script type="text/javascript">
  (function () {
    /*global _ml:true, window */
    _ml = window._ml || {};
    _ml.eid = '771';

    var s = document.getElementsByTagName('script')[0], cd = new Date(), mltag = document.createElement('script');
    mltag.type = 'text/javascript'; mltag.async = true;
    mltag.src = '//ml314.com/tag.aspx?' + cd.getDate() + cd.getMonth() + cd.getFullYear();
    s.parentNode.insertBefore(mltag, s);
  })();
</script>
<!-- End CCM Tag -->

        
        <script type="text/javascript" src="//a.fsdn.com/con/js/adframe.js?1504810730"></script>
        <script type="text/javascript">

            /*jshint ignore:start*/
            (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function() {
            (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
            m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
            })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
            /*jshint ignore:end*/

            //var $ = jQuery.noConflict(); // jshint ignore:line
        </script>
        <script type="text/javascript">
            SF.devicePixelRatio = Math.round(window.getDevicePixelRatio()*10)/10;
                ga('create', 'UA-36130941-1', {
                    'name': '1', 'sampleRate': 9});
                ga('create', 'UA-36130941-1', 'auto', {   'sampleRate': 9});
            if (typeof _ml !== 'undefined' && _ml.us) {
                if (_ml.us.tp && _ml.us.tp.length > 0) {
                    ga('set', 'dimension2', _ml.us.tp[0]);
                }
                if (_ml.us.pc && _ml.us.pc.length > 0) {
                    ga('set', 'dimension7', _ml.us.pc[0]);
                }
                ga('set', 'dimension3', _ml.us.ind);
                ga('set', 'dimension4', _ml.us.cr);
                ga('set', 'dimension5', _ml.us.cs);
                ga('set', 'dimension6', _ml.us.dm);
                ga('set', 'dimension8', _ml.us.sn);
            }

            
                ga('set', 'dimension9', 'fusioncatcher');
                ga('set', 'dimension10', 'pg_files');
                    
                ga('set', 'dimension13', 'Logged Out');
                ga('set', 'dimension14', 'No');  
                ga('set', 'dimension15', 'desktop');
                ga('set', 'dimension16', 'sync');
                ga('set', 'dimension17', SF.devicePixelRatio);
                ga('send', 'pageview');

            
        </script>
        <script type="text/javascript">
        var _paq = _paq || [];
        _paq.push(['trackPageView', document.title, {
                dimension1: 'fusioncatcher',
            dimension2: 'pg_files',
            dimension3: SF.devicePixelRatio
        }]);
        _paq.push(['enableLinkTracking']);
        (function() {
            var u="//analytics.slashdotmedia.com/";
            _paq.push(['setTrackerUrl', u+'sf.php']);
            _paq.push(['setSiteId', 39]);
            var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0];
            g.type='text/javascript'; g.async=true; g.defer=true; g.src=u+'sf.js'; s.parentNode.insertBefore(g,s);
        })();
        </script>

        <script type="text/javascript"> try{(function(){ var cb = new Date().getTime(); var s = document.createElement("script"); s.defer = true; s.src = "//tag.crsspxl.com/s1.js?d=2396&cb="+cb; var s0 = document.getElementsByTagName('script')[0]; s0.parentNode.insertBefore(s, s0); })();}catch(e){} </script>


    </head>
    <body id="pg_files" class="
    user
     bluesteel
 anonymous ">
        
        <div id="busy-spinner"></div>
        

<header id="site-header">
    <div class="wrapper">
        <a href="/" class="logo">
            <span>SourceForge</span>
        </a>
        <form method="get" action="/directory/">
            <input type="text" id="words" name="q" placeholder="Search">
        </form>
        <!--Switch to {language}-->
        <nav id="nav-site">
            <a href="/directory/" title="Browse our software.">Browse</a>
            <a href="/directory/enterprise" title="Browse our Enterprise software.">Enterprise</a>
            <a href="/blog/" title="Read the latest news from the SF HQ.">Blog</a>
            <a href="/articles/" title="Read the latest industry news about products and updates from leading cloud, network, and developer tool service providers">Articles</a>
            <a href="//deals.sourceforge.net/?utm_source=sourceforge&amp;utm_medium=navbar&amp;utm_campaign=homepage" title="Discover and Save on the Best Gear, Gadgets, and Software" class="featured-link" target="_blank">Deals</a>
            <a href="/support" title="Contact us for help and feedback.">Help</a>
            <a href="/create"  class="featured-link blue" title="Create and publish Open Source software for free.">Create</a>
        </nav>
        <nav id="nav-account">
            
              <div class="logged_out">
                  <a href="https://sourceforge.net/auth/">Log In</a>
                <span>or</span>
                <a href="https://sourceforge.net/user/registration">Join</a>
              </div>
            
        </nav>
    </div>
</header>
<header id="site-sec-header">
    <div class="wrapper">
        <nav id="nav-hubs">
            <h4>Solution Centers</h4>
            
        </nav>
        <nav id="nav-collateral">
            <a href="https://library.slashdotmedia.com/">Resources</a>
            <a href="/user/newsletters?source=sfnet_header">Newsletters</a>

            <a href="/cloud-storage-providers/?source=sfnet_header">Cloud Storage Providers</a>
            <a href="/business-voip/?source=sfnet_header">Business VoIP Providers</a>
            
            <a href="/speedtest/?source=sfnet_header">Internet Speed Test</a>
            
            <a href="/call-center-providers/?source=sfnet_header">Call Center Providers</a>
        </nav>
    </div>
</header>



        
        <header id="page-header">
    
    <div id="banner-sterling" class="sterling">
        
        
        


 
    


<div class="sticky-anchor"></div>
<div class="sticky">

<div id="SF_ProjectFiles_728x90_A_wrapped" data-id="div-gpt-ad-1393435113147-0" class="draper leaderboard  visibility_rules
 v_728_leaderboard  v_970_billboard "
>
<script type="text/javascript">
/*global googletag */
if (SF.initial_breakpoints_visible.leaderboard) {
(function(){
    
    document.write('<div id="div-gpt-ad-1393435113147-0"></div>');
}());

gptadRenderers['SF_ProjectFiles_728x90_A'] = function(){  // jshint ignore:line
    googletag.cmd.push(function() {
        googletag.display('div-gpt-ad-1393435113147-0');
    });
};
gptadRenderers['SF_ProjectFiles_728x90_A']();  // jshint ignore:line
}
</script>
</div>


</div>
<script type="text/javascript">
    if (!SF.adblock) {
        SF.Ads.stickyLeader = new SF.Stickify($('.sticky').eq(0));
    }
</script>
        
        


 
    



<div id="SF_ProjectFiles_HubIcon_200x90_A_wrapped" data-id="div-gpt-ad-1392148098424-0" class="draper hub  "
>
<script type="text/javascript">
/*global googletag */
if (true) {
(function(){
    
    document.write('<div id="div-gpt-ad-1392148098424-0"></div>');
}());

gptadRenderers['SF_ProjectFiles_HubIcon_200x90_A'] = function(){  // jshint ignore:line
    googletag.cmd.push(function() {
        googletag.display('div-gpt-ad-1392148098424-0');
    });
};
gptadRenderers['SF_ProjectFiles_HubIcon_200x90_A']();  // jshint ignore:line
}
</script>
</div>


    </div>

    
    <nav id="breadcrumbs" class="breadcrumbs">
        <ul>
            <li itemscope itemtype="http://data-vocabulary.org/Breadcrumb">
            <a itemprop="url" href="/"><span itemprop="title">Home</span></a></li>
            <li itemscope itemtype="http://data-vocabulary.org/Breadcrumb"><a itemprop="url" href="/directory"><span itemprop="title">Browse</span></a></li><li itemscope itemtype="http://data-vocabulary.org/Breadcrumb"><a itemprop="url" href="/directory/science-engineering/"><span itemprop="title">Science &amp; Engineering</span></a></li><li itemscope itemtype="http://data-vocabulary.org/Breadcrumb"><a itemprop="url" href="/directory/science-engineering/bioinformatics/"><span itemprop="title">Bio-Informatics</span></a></li><li class="project"><a href="/projects/fusioncatcher/">FusionCatcher</a></li>
            
              <li>Files</li>
            
        </ul>
    </nav>

    <div id="project-header" itemscope itemtype="http://schema.org/SoftwareApplication">
        

<section id="project-icon" class="noneditable">
    <img itemscope itemtype="http://schema.org/ImageObject" itemprop="image" alt="FusionCatcher Icon" 
    
    src="//a.fsdn.com/con/img/project_default.png"
    
 height="48" width="48"/>
</section>
<div id="project-title" class="project-mirror">
    <div>
        <div>
            <h1 itemprop="name">
                    <a href="/projects/fusioncatcher/" itemprop="url">FusionCatcher</a></h1>
                    <span id="dev-status" class="beta">beta</span>
        </div>
        <h2>
            Somatic fusion-genes finder for RNA-seq data
        </h2>
        <p id="maintainers" itemprop="author" itemscope itemtype="http://schema.org/Person">
        Brought to you by:
        <a href="/u/daninico/" itemprop="url"><span itemprop="name">daninico</span></a></p>
        
    </div>
    

</div>
<nav id="project-nav">
    
    <div id="top_nav"><div id="top_nav_admin">
        <ul class="dropdown">
            
            <li >
                <a href="/projects/fusioncatcher/?source=navbar"
                >
                <span>Summary</span></a>
                
            </li>
            
            <li class="selected">
                <a href="/projects/fusioncatcher/files/?source=navbar"
                >
                <span>Files</span></a>
                
            </li>
            
            <li >
                <a href="/projects/fusioncatcher/reviews?source=navbar"
                >
                <span>Reviews</span></a>
                
            </li>
            
            <li >
                <a href="/projects/fusioncatcher/support?source=navbar"
                >
                <span>Support</span></a>
                
            </li>
            
            <li >
                <a href="/p/fusioncatcher/wiki/?source=navbar"
                >
                <span>Wiki</span></a>
                
            </li>
            
            <li >
                <a href="/p/fusioncatcher/code/?source=navbar"
                >
                <span>Code</span></a>
                
            </li>
            
            <li >
                <a href="/p/fusioncatcher/tickets/?source=navbar"
                >
                <span>Tickets</span></a>
                
            </li>
            
            <li >
                <a href="/p/fusioncatcher/discussion/?source=navbar"
                >
                <span>Discussion</span></a>
                
            </li>
            
            
        </ul>
        
    </div></div>
    
</nav>


    </div>
</header>


        
<div id="messages">
</div>



        <div id="page-body">
    <noscript>
        <p>The interactive file manager requires Javascript. Please enable it or use <a href="https://sourceforge.net/p/forge/documentation/Release%20Files%20for%20Download#scp">sftp or scp</a>.
        <br/>You may still <em>browse</em> the files here.</p>
    </noscript>
    <div id="files">
      <div class="download-bar">Looking for the latest version? <strong>
            <a href="/projects/fusioncatcher/files/latest/download?source=files" title="/extra/setuptools-36.2.7.zip:  released on 2017-08-07 05:16:55 UTC">
                <span>Download setuptools-36.2.7.zip (716.4 kB)</span>
            </a>
            </strong>
        </div><div class="files-breadcrumb">
            <span class="actions"><a href="/projects/fusioncatcher/rss?path=/" title="Monitor this project's files using RSS" data-icon="f" class="ico ico-feed"></a></span>
            
            Home
            
            

        </div>

        <table id="files_list" title="This is a list of the files for this project">
            <col class="name-column">
            <col class="date-column">
            <col class="size-column">
            <col class="downloads-column">
            <col class="status-column">
            <thead>
                <tr>
                    <th title="The file or folder's name" id="files_name_h" class="first">Name</th>
                    <th title="The file or folder's last modified date" id="files_date_h" class="opt">Modified</th>
                    <th title="The file size" id="files_size_h" class="opt">Size</th>
                    <th title="The weekly download count" id="files_downloads_h" class="opt">Downloads / Week</th>
                    <th title="A snapshot of the file or folder's status" id="files_status_h" class="typesort">Status</th>
                </tr>
                
            </thead>
            <tbody>
                
                <tr title="extra" class="folder ">
                    <th scope="row" headers="files_name_h"><span data-icon="o" class="ico ico-folder"></span>
                        
                        
                        <a href="/projects/fusioncatcher/files/extra/"
                           title="Click to enter extra"
                           class="name">
                        extra</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-08-07 05:16:55 UTC">2017-08-07</abbr></td>
                    <td headers="files_size_h" class="opt"></td>
                    <td headers="files_downloads_h" class="opt">
                        2<a href="/projects/fusioncatcher/files/extra/stats/timeline" class="fs-stats ui-corner-all folder" rel="nofollow"><img class="fs-sparklines folder" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:0,0,2,0,0,0,0&amp;chds=0,36" title="2 weekly downloads" alt="2 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status folder"></td>
                </tr>
                <tr title="test" class="folder ">
                    <th scope="row" headers="files_name_h"><span data-icon="o" class="ico ico-folder"></span>
                        
                        
                        <a href="/projects/fusioncatcher/files/test/"
                           title="Click to enter test"
                           class="name">
                        test</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-06-29 09:35:09 UTC">2017-06-29</abbr></td>
                    <td headers="files_size_h" class="opt"></td>
                    <td headers="files_downloads_h" class="opt">
                        8<a href="/projects/fusioncatcher/files/test/stats/timeline" class="fs-stats ui-corner-all folder" rel="nofollow"><img class="fs-sparklines folder" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:2,1,0,5,0,0,0&amp;chds=0,36" title="8 weekly downloads" alt="8 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status folder"></td>
                </tr>
                <tr title="data" class="folder ">
                    <th scope="row" headers="files_name_h"><span data-icon="o" class="ico ico-folder"></span>
                        
                        
                        <a href="/projects/fusioncatcher/files/data/"
                           title="Click to enter data"
                           class="name">
                        data</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-06-26 10:51:10 UTC">2017-06-26</abbr></td>
                    <td headers="files_size_h" class="opt"></td>
                    <td headers="files_downloads_h" class="opt">
                        79<a href="/projects/fusioncatcher/files/data/stats/timeline" class="fs-stats ui-corner-all folder" rel="nofollow"><img class="fs-sparklines folder" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:36,27,11,5,0,0,0&amp;chds=0,36" title="79 weekly downloads" alt="79 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status folder"></td>
                </tr>
                <tr title="old" class="folder ">
                    <th scope="row" headers="files_name_h"><span data-icon="o" class="ico ico-folder"></span>
                        
                        
                        <a href="/projects/fusioncatcher/files/old/"
                           title="Click to enter old"
                           class="name">
                        old</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-06-22 10:12:19 UTC">2017-06-22</abbr></td>
                    <td headers="files_size_h" class="opt"></td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads folder" href="/projects/fusioncatcher/files/old/stats/timeline" title="8 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status folder"></td>
                </tr>
                <tr title="examples" class="folder ">
                    <th scope="row" headers="files_name_h"><span data-icon="o" class="ico ico-folder"></span>
                        
                        
                        <a href="/projects/fusioncatcher/files/examples/"
                           title="Click to enter examples"
                           class="name">
                        examples</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-02-17 07:31:07 UTC">2017-02-17</abbr></td>
                    <td headers="files_size_h" class="opt"></td>
                    <td headers="files_downloads_h" class="opt">
                        5<a href="/projects/fusioncatcher/files/examples/stats/timeline" class="fs-stats ui-corner-all folder" rel="nofollow"><img class="fs-sparklines folder" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:0,4,0,1,0,0,0&amp;chds=0,36" title="5 weekly downloads" alt="5 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status folder"></td>
                </tr>
                <tr title="bootstrap.py" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/bootstrap.py/download"
                           title="Click to download bootstrap.py"
                           class="name">
                        bootstrap.py</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-08-07 05:31:49 UTC">2017-08-07</abbr></td>
                    <td headers="files_size_h" class="opt">96.5 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        19<a href="/projects/fusioncatcher/files/bootstrap.py/stats/timeline" class="fs-stats ui-corner-all file" rel="nofollow"><img class="fs-sparklines file" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:7,6,3,2,1,0,0&amp;chds=0,36" title="19 weekly downloads" alt="19 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.7c.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7c.zip/download"
                           title="Click to download fusioncatcher_v0.99.7c.zip"
                           class="name">
                        fusioncatcher_v0.99.7c.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-06-30 07:15:27 UTC">2017-06-30</abbr></td>
                    <td headers="files_size_h" class="opt">2.7 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        25<a href="/projects/fusioncatcher/files/fusioncatcher_v0.99.7c.zip/stats/timeline" class="fs-stats ui-corner-all file" rel="nofollow"><img class="fs-sparklines file" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:11,7,4,2,1,0,0&amp;chds=0,36" title="25 weekly downloads" alt="25 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="README.md" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/README.md/download"
                           title="Click to download README.md"
                           class="name">
                        README.md</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-06-29 12:27:59 UTC">2017-06-29</abbr></td>
                    <td headers="files_size_h" class="opt">3.3 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        1<a href="/projects/fusioncatcher/files/README.md/stats/timeline" class="fs-stats ui-corner-all file" rel="nofollow"><img class="fs-sparklines file" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:0,0,0,1,0,0,0&amp;chds=0,36" title="1 weekly downloads" alt="1 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.7b.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7b.zip/download"
                           title="Click to download fusioncatcher_v0.99.7b.zip"
                           class="name">
                        fusioncatcher_v0.99.7b.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-03-17 08:12:00 UTC">2017-03-17</abbr></td>
                    <td headers="files_size_h" class="opt">2.6 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.7b.zip/stats/timeline" title="683 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.7a.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7a.zip/download"
                           title="Click to download fusioncatcher_v0.99.7a.zip"
                           class="name">
                        fusioncatcher_v0.99.7a.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2017-02-16 07:47:36 UTC">2017-02-16</abbr></td>
                    <td headers="files_size_h" class="opt">2.6 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.7a.zip/stats/timeline" title="61 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.6a.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip/download"
                           title="Click to download fusioncatcher_v0.99.6a.zip"
                           class="name">
                        fusioncatcher_v0.99.6a.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2016-05-20 11:57:50 UTC">2016-05-20</abbr></td>
                    <td headers="files_size_h" class="opt">2.6 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip/stats/timeline" title="1,730 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.5a.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.5a.zip/download"
                           title="Click to download fusioncatcher_v0.99.5a.zip"
                           class="name">
                        fusioncatcher_v0.99.5a.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2016-01-20 10:42:14 UTC">2016-01-20</abbr></td>
                    <td headers="files_size_h" class="opt">2.5 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.5a.zip/stats/timeline" title="571 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.4e.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4e.zip/download"
                           title="Click to download fusioncatcher_v0.99.4e.zip"
                           class="name">
                        fusioncatcher_v0.99.4e.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2015-12-09 10:14:42 UTC">2015-12-09</abbr></td>
                    <td headers="files_size_h" class="opt">2.5 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.4e.zip/stats/timeline" title="231 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.4d.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4d.zip/download"
                           title="Click to download fusioncatcher_v0.99.4d.zip"
                           class="name">
                        fusioncatcher_v0.99.4d.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2015-10-19 10:52:24 UTC">2015-10-19</abbr></td>
                    <td headers="files_size_h" class="opt">2.4 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.4d.zip/stats/timeline" title="500 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.4c.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4c.zip/download"
                           title="Click to download fusioncatcher_v0.99.4c.zip"
                           class="name">
                        fusioncatcher_v0.99.4c.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2015-06-26 05:35:08 UTC">2015-06-26</abbr></td>
                    <td headers="files_size_h" class="opt">2.4 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.4c.zip/stats/timeline" title="386 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.4b.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4b.zip/download"
                           title="Click to download fusioncatcher_v0.99.4b.zip"
                           class="name">
                        fusioncatcher_v0.99.4b.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2015-05-15 09:25:30 UTC">2015-05-15</abbr></td>
                    <td headers="files_size_h" class="opt">2.4 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.4b.zip/stats/timeline" title="207 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.4a.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4a.zip/download"
                           title="Click to download fusioncatcher_v0.99.4a.zip"
                           class="name">
                        fusioncatcher_v0.99.4a.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2015-03-18 06:35:25 UTC">2015-03-18</abbr></td>
                    <td headers="files_size_h" class="opt">2.3 MB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.4a.zip/stats/timeline" title="248 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.3e.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3e.zip/download"
                           title="Click to download fusioncatcher_v0.99.3e.zip"
                           class="name">
                        fusioncatcher_v0.99.3e.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-12-11 07:13:16 UTC">2014-12-11</abbr></td>
                    <td headers="files_size_h" class="opt">650.4 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.3e.zip/stats/timeline" title="535 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.3d.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3d.zip/download"
                           title="Click to download fusioncatcher_v0.99.3d.zip"
                           class="name">
                        fusioncatcher_v0.99.3d.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-11-08 08:27:29 UTC">2014-11-08</abbr></td>
                    <td headers="files_size_h" class="opt">644.7 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.3d.zip/stats/timeline" title="68 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.3c.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3c.zip/download"
                           title="Click to download fusioncatcher_v0.99.3c.zip"
                           class="name">
                        fusioncatcher_v0.99.3c.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-10-20 06:44:40 UTC">2014-10-20</abbr></td>
                    <td headers="files_size_h" class="opt">623.4 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.3c.zip/stats/timeline" title="173 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.3b.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3b.zip/download"
                           title="Click to download fusioncatcher_v0.99.3b.zip"
                           class="name">
                        fusioncatcher_v0.99.3b.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-07-29 10:59:05 UTC">2014-07-29</abbr></td>
                    <td headers="files_size_h" class="opt">499.8 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.3b.zip/stats/timeline" title="390 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.3a.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3a.zip/download"
                           title="Click to download fusioncatcher_v0.99.3a.zip"
                           class="name">
                        fusioncatcher_v0.99.3a.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-05-15 05:25:57 UTC">2014-05-15</abbr></td>
                    <td headers="files_size_h" class="opt">428.6 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.3a.zip/stats/timeline" title="215 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.2.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.2.zip/download"
                           title="Click to download fusioncatcher_v0.99.2.zip"
                           class="name">
                        fusioncatcher_v0.99.2.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-07 07:14:37 UTC">2014-02-07</abbr></td>
                    <td headers="files_size_h" class="opt">414.4 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        1<a href="/projects/fusioncatcher/files/fusioncatcher_v0.99.2.zip/stats/timeline" class="fs-stats ui-corner-all file" rel="nofollow"><img class="fs-sparklines file" src="https://a.fsdn.com/con/img/gallery/opa75.png" data-src="https://chart.googleapis.com/chart?cht=ls&amp;chs=20x16&amp;chco=0685c6&amp;chm=B,0685c6,0,0,0&amp;chd=t:0,0,0,0,1,0,0&amp;chds=0,36" title="1 weekly downloads" alt="1 weekly downloads" width="20" height="16" /></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.1.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.1.zip/download"
                           title="Click to download fusioncatcher_v0.99.1.zip"
                           class="name">
                        fusioncatcher_v0.99.1.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 08:03:30 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">379.7 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.1.zip/stats/timeline" title="2 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.99.0.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.0.zip/download"
                           title="Click to download fusioncatcher_v0.99.0.zip"
                           class="name">
                        fusioncatcher_v0.99.0.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 07:29:35 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">367.9 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.99.0.zip/stats/timeline" title="1 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.95.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.95.zip/download"
                           title="Click to download fusioncatcher_v0.95.zip"
                           class="name">
                        fusioncatcher_v0.95.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 06:49:44 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">239.9 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.95.zip/stats/timeline" title="1 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.96.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.96.zip/download"
                           title="Click to download fusioncatcher_v0.96.zip"
                           class="name">
                        fusioncatcher_v0.96.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 06:49:44 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">287.4 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.96.zip/stats/timeline" title="2 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.93.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.93.zip/download"
                           title="Click to download fusioncatcher_v0.93.zip"
                           class="name">
                        fusioncatcher_v0.93.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 06:49:43 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">206.8 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.93.zip/stats/timeline" title="1 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.94.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.94.zip/download"
                           title="Click to download fusioncatcher_v0.94.zip"
                           class="name">
                        fusioncatcher_v0.94.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2014-02-01 06:49:43 UTC">2014-02-01</abbr></td>
                    <td headers="files_size_h" class="opt">247.4 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.94.zip/stats/timeline" title="3 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.98.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.98.zip/download"
                           title="Click to download fusioncatcher_v0.98.zip"
                           class="name">
                        fusioncatcher_v0.98.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2013-06-27 06:02:45 UTC">2013-06-27</abbr></td>
                    <td headers="files_size_h" class="opt">362.0 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.98.zip/stats/timeline" title="69 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
                <tr title="fusioncatcher_v0.97.zip" class="file ">
                    <th scope="row" headers="files_name_h">
                        
                        
                        <a href="https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.97.zip/download"
                           title="Click to download fusioncatcher_v0.97.zip"
                           class="name">
                        fusioncatcher_v0.97.zip</a></th>
                    <td headers="files_date_h" class="opt"><abbr title="2013-06-24 18:03:11 UTC">2013-06-24</abbr></td>
                    <td headers="files_size_h" class="opt">304.3 kB</td>
                    <td headers="files_downloads_h" class="opt">
                        0<a class="fs-stats ui-corner-all no-recent-downloads file" href="/projects/fusioncatcher/files/fusioncatcher_v0.97.zip/stats/timeline" title="4 downloads (all-time), none recently." rel="nofollow"></a></td>
                    <td headers="files_status_h" class="status file"></td>
                </tr>
            </tbody>
            <tfoot>
                <tr>
                    <td id="totals"><span class="label">Totals: </span>31 Items</td>
                    <td class="opt">&nbsp;</td>
                    <td headers="files_size_h" class="opt">30.8 MB</td>
                    <td headers="files_downloads_h" class="opt">140</td>
                    <td headers="files_status_h">
                        
                    </td>
                </tr>
            </tfoot>
        </table>
        <div id="files-drawer" class="fs-widget fs-drawer consumer">
        </div>
        <div id="readme">
            
            <div class="content format-markdown"><h1>FusionCatcher</h1>
<p>Finder of somatic fusion-genes in RNA-seq data.</p>
<h2>Download / Install / Update / Upgrade <a href="http://github.com/ndaniel/fusioncatcher">FusionCatcher</a></h2>
<p>Use this one-line command:</p>
<p><code>wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py &amp;&amp; python bootstrap.py -t --download</code></p>
<p>If one wants to have all the questions asked by boostrap.py answered automatically with yes then add <code>-y</code> to the 
command above. For more installing options, see:</p>
<p><code>bootstrap.py --help</code></p>
<p>On Ubuntu Linux running this command before installing FusionCatcher using <code>bootstrap.py</code> would help making the installation process smoother:</p>
<p><code>sudo apt-get install wget gawk gcc g++ make cmake automake curl unzip zip bzip2 tar gzip pigz parallel build-essential libncurses5-dev libc6-dev zlib1g zlib1g-dev libtbb-dev libtbb2 python python-dev python-numpy python-biopython python-xlrd python-openpyxl default-jdk</code></p>
<h2>Description</h2>
<p>FusionCatcher searches for novel/known somatic fusion genes, translocations, and
chimeras in RNA-seq data (paired-end or single-end reads from Illumina NGS platforms 
like Solexa/HiSeq/NextSeq/MiSeq/MiniSeq) from diseased samples.</p>
<p>The aims of FusionCatcher are:
 * very good detection rate for finding candidate somatic fusion
   genes (see somatic mutations; using a matched normal sample is
   optional; several databases of known fusion genes found in healthy
   samples are used as a list of known false positives; biological
   knowledge is used, like for example gene fusion between a gene and
   its pseudogene is filtered out),
 * very good RT-PCR validation rate of found candidate somatic fusion
   genes (this is very important for us),
 * very easy to use (i.e. no a priori knowledge of bioinformatic
   databases and bioinformatics is needed in order to run FusionCatcher BUT
   Linux/Unix knowledge is needed; it allows a very high level of control
   for expert users),
 * to be as automatic as possible (i.e. the FusionCatcher will choose
   automatically the best parameters in order to find candidate somatic
   fusion genes, e.g. finding automatically the adapters, quality trimming
   of reads, building the exon-exon junctions automatically based on the
   length of the reads given as input, etc. while giving also full control
   to expert users) while providing the best possible detection rate for
   finding somatic fusion genes (with a very low rate of false positives
   but a very good precision).</p>
<h2>Manual</h2>
<p>A detailed manual is available <a href="doc/manual.md">here</a>.</p>
<h2>Forum</h2>
<p>A forum for FusionCatcher is available at 
<a href="http://groups.google.com/d/forum/fusioncatcher">Google Groups</a>.</p>
<h2>Release history</h2>
<p>A complete release history can be found <a href="NEWS">here</a>.</p>
<h2>Official releases</h2>
<p>Old releases and the latest official release of FusionCatcher are on <a href="https://sourceforge.net/projects/fusioncatcher/files/">https://sourceforge.net/projects/fusioncatcher/files/</a></p>
<h2>Citing</h2>
<p>D. Nicorici, M. Satalan, H. Edgren, S. Kangaspeska, A. Murumagi, O. Kallioniemi,
S. Virtanen, O. Kilkku, FusionCatcher – a tool for finding somatic fusion genes
in paired-end RNA-sequencing data, bioRxiv, Nov. 2014, 
<a href="http://dx.doi.org/10.1101/011650">DOI:10.1101/011650</a></p></div>
            <span class="meta">Source: README.md, updated 2017-06-29</span>
            
        </div>
    </div>

    
    <div id="files-sidebar" class="scroll-fixable" data-floor-compensate="145">
        <div class="sterling">
        


 
    
    



<div id="SF_ProjectFiles_300x250_A_wrapped" data-id="div-gpt-ad-1392147725721-0" class="draper multisize  "
>
<script type="text/javascript">
/*global googletag */
if (true) {
(function(){
    
    document.write('<div id="div-gpt-ad-1392147725721-0"></div>');
}());

gptadRenderers['SF_ProjectFiles_300x250_A'] = function(){  // jshint ignore:line
    googletag.cmd.push(function() {
        googletag.display('div-gpt-ad-1392147725721-0');
    });
};
gptadRenderers['SF_ProjectFiles_300x250_A']();  // jshint ignore:line
}
</script>
</div>


        </div>
                <aside class="sidebar-widget">
        <header>
        <h3>Recommended Projects</h3>
        </header>
    <ul class="vertical">

        

            
                <li class="item odd">
                <a href="/projects/bowtie-bio/?source=recommended" title="Bowtie"><img alt="Icon" src="//a.fsdn.com/con/img/project_default.png" /></a>
                <div class="pinfo-content recommended">
                <a class="project-name" href="/projects/bowtie-bio/?source=recommended"
                        title="Learn more about Bowtie ">Bowtie</a>
                        <div class="summary">
                        
                        </div>
                    </div>

            </li>
            
                <li class="item even">
                <a href="/projects/rseqc/?source=recommended" title="RSeQC"><img alt="Icon" src="//a.fsdn.com/con/img/project_default.png" /></a>
                <div class="pinfo-content recommended">
                <a class="project-name" href="/projects/rseqc/?source=recommended"
                        title="Learn more about RSeQC ">RSeQC</a>
                        <div class="summary">
                        RNA-seq data QC
                        </div>
                    </div>

            </li>
            
                <li class="item odd last">
                <a href="/projects/samtools/?source=recommended" title="SAM tools"><img alt="Icon" src="//a.fsdn.com/con/img/project_default.png" /></a>
                <div class="pinfo-content recommended">
                <a class="project-name" href="/projects/samtools/?source=recommended"
                        title="Learn more about SAM tools ">SAM tools</a>
                        <div class="summary">
                        
                        </div>
                    </div>

            </li>
            
        
    </ul>
</aside>


        


 
    



<div id="SF_ProjectFiles_300x250_B_wrapped" data-id="div-gpt-ad-1392148208789-0" class="draper medrec  "
>
<script type="text/javascript">
/*global googletag */
if (true) {
(function(){
    
    document.write('<div id="div-gpt-ad-1392148208789-0"></div>');
}());

gptadRenderers['SF_ProjectFiles_300x250_B'] = function(){  // jshint ignore:line
    googletag.cmd.push(function() {
        googletag.display('div-gpt-ad-1392148208789-0');
    });
};
gptadRenderers['SF_ProjectFiles_300x250_B']();  // jshint ignore:line
}
</script>
</div>



        
            <div class="sterling" id="deals-widget">
                
    <aside class="sidebar-widget stackcommerce-widget">
            <header id="stackcommerce-header">
                <img src="//a.fsdn.com/con/img/sf_logo250x32.png" width="125" height="16" alt="sourceforge.net logo" />
                <h3>Deals</h3>
            </header>
            


 
    



<div id="SF-300x250_wrapped" data-id="div-gpt-ad-1392148208795-0" class="draper medrec  visibility_rules
 v_300_large "
>
<script type="text/javascript">
/*global googletag */
if (SF.initial_breakpoints_visible.large) {
(function(){
    
    document.write('<div id="div-gpt-ad-1392148208795-0"></div>');
}());

gptadRenderers['SF-300x250'] = function(){  // jshint ignore:line
    googletag.cmd.push(function() {
        googletag.display('div-gpt-ad-1392148208795-0');
    });
};
gptadRenderers['SF-300x250']();  // jshint ignore:line
}
</script>
</div>


    </aside>
    

            </div>
        
    </div>
    
    
<script type="text/javascript">
if (!SF.adblock) {
    
    SF.Ads.scrollFixableEnabled = true;
    
}
</script>

        </div>
        

        
            
    <div id="overlay-blockthis-wrapper">
    <div id="overlay-blockthis">
        <div class="instructions">
            <h2>Thanks for helping keep SourceForge clean.</h2>
            <p>
            <u>Screenshot instructions:</u><br>
            <a data-external target=_blank href="http://windows.microsoft.com/en-us/windows/take-screen-capture-print-screen#take-screen-capture-print-screen=windows-8">Windows</a><br>
            <a data-external target=_blank href="https://support.apple.com/en-us/HT201361">Mac</a><br>
            <a data-external target=_blank href="https://access.redhat.com/solutions/2178">Red Hat Linux</a> &nbsp;
            <a data-external target=_blank href="https://help.ubuntu.com/stable/ubuntu-help/screen-shot-record.html">Ubuntu</a>
            </p>
            <p>
                <u>Click URL instructions:</u><br>
                Right-click on ad, choose "Copy Link", then paste here &rarr;<br>
                (This may not be possible with some types of ads)
            </p>
            <a class="more-info" href="https://sourceforge.net/p/forge/documentation/Report%20a%20problem%20with%20Ad%20content/" target="_blank">More information about our ad policies</a>
        </div>
        <form class="dropzone" action="/api/instrumentation/gpt" id="blockthisForm" method="POST">
            <a href="#" id="btn-blockthis-close">X</a>
            
  <input type="hidden" name="_visit_cookie" value=""/>
                <input type='hidden' name='timestamp' value='1505119436'/>
                <input type='hidden' name='spinner' value='XbzIVz35oGfVNAzctHsU22_rkcfw'/>
                <p class='H9effeddc4f6fef7a225a171e663950d2074984de'><label for='XaFp6oRsRKcXVhaERifGdrF_pnlw'>You seem to have CSS turned off.
             Please don't fill out this field.</label><input id='XaFp6oRsRKcXVhaERifGdrF_pnlw' name='XaVp6oRsRKW3LlQu6Km5BfvcL0VM' type=
             'text'/></p>
                <p class='H9effeddc4f6fef7a225a171e663950d2074984de'><label for='XaFp6oRsRKcTVhaERifGdrF_pnlw'>You seem to have CSS turned off.
             Please don't fill out this field.</label><input id='XaFp6oRsRKcTVhaERifGdrF_pnlw' name='XaVp6oRsRKG3LlQu6Km5BfvcL0VM' type=
             'text'/></p>
            <p>Briefly describe the problem (required):
            <input name="XZFZwvB0acIU5alhDhkOg523Q2os" type="text" required>
            </p>

            <div class="upload-text">Upload screenshot of ad (required):</div>
            <div id='upload-it'>
                <a href="#" onclick="return false">Select a file</a>, or drag & drop file here.
            </div>
            <div id="upload-it-placeholder"></div> 

            <div class="dropzone-previews" style="display: none"></div>
            <div class="dz-message" style="display: none"></div> 
            
            <div id="dropzone-preview-template" style="display: none">
                <div class="dz-preview dz-file-preview">
                    <img data-dz-thumbnail src="data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs=" alt=""/>
                    <div class="dz-success-mark"><span>✔</span></div>
                    <div class="dz-error-mark"><span>✘</span></div>
                    <div class="dz-error-message"><span data-dz-errormessage></span></div>
                </div>
            </div>
            <p></p>
            <p>Please provide the ad click URL, if possible:
            <input name="XZlF5ph0DRoA_b6:riPmh71GT1PE" type="url">
            </p>
            <textarea id="gpt-info" name="Xa1Z0ux_wn2NxlAOGaWA7NFpL4pU"></textarea>
            <input type="submit" id="btn-blockthis-submit" value="Submit Report">
        </form>
    </div>
    </div>

        

        <script type="text/javascript" src="//ads.pro-market.net/ads/scripts/site-143572.js"></script>

<footer id="site-footer">
    <div class="wrapper">
        <nav>
            <h5>SourceForge</h5>
            <a href="/about">About</a>
            <a href="/blog/category/sitestatus/">Site Status</a>
            <a href="http://twitter.com/sfnet_ops">@sfnet_ops</a>
        </nav>
        <nav>
            <h5>Find and Develop Software</h5>
            <a href="/create/">Create a Project</a>
            <a href="/directory/">Software Directory</a>
            <a href="/top">Top Downloaded Projects</a>
        </nav>
        <nav>
            <h5>Community</h5>
            <a href="/blog/">Blog</a>
            <a href="http://twitter.com/sourceforge">@sourceforge</a>
            <a href="https://library.slashdotmedia.com/">Resources</a>
        </nav>
        <nav>
            <h5>Help</h5>
            <a href="http://p.sf.net/sourceforge/docs">Site Documentation</a>
            <a href="/support">Support Request</a>
        </nav>
    </div>
</footer>
<footer id="site-copyright-footer">
    <div class="wrapper">
        <div id="copyright">
            &copy; 2017 Slashdot Media. All Rights Reserved.<br />
        </div>
        <nav>
            <a href="https://slashdotmedia.com/terms-of-use/">Terms</a>
            <a href="https://slashdotmedia.com/privacy-statement/">Privacy</a>
            <span id='teconsent'></span>
            <a href="https://slashdotmedia.com/opt-out-choices/">Opt Out Choices</a>
            <a href="https://slashdotmedia.com">Advertise</a>
        </nav>
    </div>
</footer>

        
<div id="newsletter-floating" class="goth-form">
    <h2>Get latest updates about Open Source Projects, Conferences and News.</h2>
    <div class="teaser">
    <input type="submit" value="Sign Up" class="bt">
    </div>
    <div class="form-itself">
    


<form action="/user/newsletters/subscribe" method="post" class="newsletter-subscribe-form compliance-form "  >

    <div class="form">
        <div class="fielderror"></div>
        <input type="email" name="Xald4rhcEgXPbP6AZtbKT1hVE3m8"  placeholder="email@address.com" value="" required>

        <label class="input-set stacked input-set-country">
            <span class="label">Country</span>
            <span class="input">
<select id="country-floating" name="XaFF6uhAca4zVhaERifGdrF_pnlw" required >
    <option value=""></option>
    <option value="AF">Afghanistan</option>
    <option value="AX">Aland Islands</option>
    <option value="AL">Albania</option>
    <option value="DZ">Algeria</option>
    <option value="AS">American Samoa</option>
    <option value="AD">Andorra</option>
    <option value="AO">Angola</option>
    <option value="AI">Anguilla</option>
    <option value="AQ">Antarctica</option>
    <option value="AG">Antigua and Barbuda</option>
    <option value="AR">Argentina</option>
    <option value="AM">Armenia</option>
    <option value="AW">Aruba</option>
    <option value="AU">Australia</option>
    <option value="AT">Austria</option>
    <option value="AZ">Azerbaijan</option>
    <option value="BS">Bahamas</option>
    <option value="BH">Bahrain</option>
    <option value="BD">Bangladesh</option>
    <option value="BB">Barbados</option>
    <option value="BY">Belarus</option>
    <option value="BE">Belgium</option>
    <option value="BZ">Belize</option>
    <option value="BJ">Benin</option>
    <option value="BM">Bermuda</option>
    <option value="BT">Bhutan</option>
    <option value="BO">Bolivia</option>
    <option value="BA">Bosnia and Herzegovina</option>
    <option value="BW">Botswana</option>
    <option value="BV">Bouvet Island</option>
    <option value="BR">Brazil</option>
    <option value="IO">British Indian Ocean Territory</option>
    <option value="BN">Brunei Darussalam</option>
    <option value="BG">Bulgaria</option>
    <option value="BF">Burkina Faso</option>
    <option value="BI">Burundi</option>
    <option value="KH">Cambodia</option>
    <option value="CM">Cameroon</option>
    <option value="CA">Canada</option>
    <option value="CV">Cape Verde</option>
    <option value="KY">Cayman Islands</option>
    <option value="CF">Central African Republic</option>
    <option value="TD">Chad</option>
    <option value="CL">Chile</option>
    <option value="CN">China</option>
    <option value="CX">Christmas Island</option>
    <option value="CC">Cocos (Keeling) Islands</option>
    <option value="CO">Colombia</option>
    <option value="KM">Comoros</option>
    <option value="CG">Congo</option>
    <option value="CD">Congo, The Democratic Republic of the</option>
    <option value="CK">Cook Islands</option>
    <option value="CR">Costa Rica</option>
    <option value="CI">Cote D&#39;Ivoire</option>
    <option value="HR">Croatia</option>
    <option value="CU">Cuba</option>
    <option value="CY">Cyprus</option>
    <option value="CZ">Czech Republic</option>
    <option value="DK">Denmark</option>
    <option value="DJ">Djibouti</option>
    <option value="DM">Dominica</option>
    <option value="DO">Dominican Republic</option>
    <option value="EC">Ecuador</option>
    <option value="EG">Egypt</option>
    <option value="SV">El Salvador</option>
    <option value="GQ">Equatorial Guinea</option>
    <option value="ER">Eritrea</option>
    <option value="EE">Estonia</option>
    <option value="ET">Ethiopia</option>
    <option value="FK">Falkland Islands (Malvinas)</option>
    <option value="FO">Faroe Islands</option>
    <option value="FJ">Fiji</option>
    <option value="FI">Finland</option>
    <option value="FR">France</option>
    <option value="GF">French Guiana</option>
    <option value="PF">French Polynesia</option>
    <option value="TF">French Southern Territories</option>
    <option value="GA">Gabon</option>
    <option value="GM">Gambia</option>
    <option value="GE">Georgia</option>
    <option value="DE">Germany</option>
    <option value="GH">Ghana</option>
    <option value="GI">Gibraltar</option>
    <option value="GR">Greece</option>
    <option value="GL">Greenland</option>
    <option value="GD">Grenada</option>
    <option value="GP">Guadeloupe</option>
    <option value="GU">Guam</option>
    <option value="GT">Guatemala</option>
    <option value="GG">Guernsey</option>
    <option value="GN">Guinea</option>
    <option value="GW">Guinea-Bissau</option>
    <option value="GY">Guyana</option>
    <option value="HT">Haiti</option>
    <option value="HM">Heard Island and McDonald Islands</option>
    <option value="VA">Holy See (Vatican City State)</option>
    <option value="HN">Honduras</option>
    <option value="HK">Hong Kong</option>
    <option value="HU">Hungary</option>
    <option value="IS">Iceland</option>
    <option value="IN">India</option>
    <option value="ID">Indonesia</option>
    <option value="IR">Iran, Islamic Republic of</option>
    <option value="IQ">Iraq</option>
    <option value="IE">Ireland</option>
    <option value="IM">Isle of Man</option>
    <option value="IL">Israel</option>
    <option value="IT">Italy</option>
    <option value="JM">Jamaica</option>
    <option value="JP">Japan</option>
    <option value="JE">Jersey</option>
    <option value="JO">Jordan</option>
    <option value="KZ">Kazakhstan</option>
    <option value="KE">Kenya</option>
    <option value="KI">Kiribati</option>
    <option value="KP">Korea, Democratic People&#39;s Republic of</option>
    <option value="KR">Korea, Republic of</option>
    <option value="XK">Kosovo</option>
    <option value="KW">Kuwait</option>
    <option value="KG">Kyrgyzstan</option>
    <option value="LA">Lao People&#39;s Democratic Republic</option>
    <option value="LV">Latvia</option>
    <option value="LB">Lebanon</option>
    <option value="LS">Lesotho</option>
    <option value="LR">Liberia</option>
    <option value="LY">Libyan Arab Jamahiriya</option>
    <option value="LI">Liechtenstein</option>
    <option value="LT">Lithuania</option>
    <option value="LU">Luxembourg</option>
    <option value="MO">Macau</option>
    <option value="MK">Macedonia</option>
    <option value="MG">Madagascar</option>
    <option value="MW">Malawi</option>
    <option value="MY">Malaysia</option>
    <option value="MV">Maldives</option>
    <option value="ML">Mali</option>
    <option value="MT">Malta</option>
    <option value="MH">Marshall Islands</option>
    <option value="MQ">Martinique</option>
    <option value="MR">Mauritania</option>
    <option value="MU">Mauritius</option>
    <option value="YT">Mayotte</option>
    <option value="MX">Mexico</option>
    <option value="FM">Micronesia, Federated States of</option>
    <option value="MD">Moldova, Republic of</option>
    <option value="MC">Monaco</option>
    <option value="MN">Mongolia</option>
    <option value="ME">Montenegro</option>
    <option value="MS">Montserrat</option>
    <option value="MA">Morocco</option>
    <option value="MZ">Mozambique</option>
    <option value="MM">Myanmar</option>
    <option value="NA">Namibia</option>
    <option value="NR">Nauru</option>
    <option value="NP">Nepal</option>
    <option value="NL">Netherlands</option>
    <option value="AN">Netherlands Antilles</option>
    <option value="NC">New Caledonia</option>
    <option value="NZ">New Zealand</option>
    <option value="NI">Nicaragua</option>
    <option value="NE">Niger</option>
    <option value="NG">Nigeria</option>
    <option value="NU">Niue</option>
    <option value="NF">Norfolk Island</option>
    <option value="MP">Northern Mariana Islands</option>
    <option value="NO">Norway</option>
    <option value="OM">Oman</option>
    <option value="PK">Pakistan</option>
    <option value="PW">Palau</option>
    <option value="PS">Palestinian Territory</option>
    <option value="PA">Panama</option>
    <option value="PG">Papua New Guinea</option>
    <option value="PY">Paraguay</option>
    <option value="PE">Peru</option>
    <option value="PH">Philippines</option>
    <option value="PN">Pitcairn Islands</option>
    <option value="PL">Poland</option>
    <option value="PT">Portugal</option>
    <option value="PR">Puerto Rico</option>
    <option value="QA">Qatar</option>
    <option value="RE">Reunion</option>
    <option value="RO">Romania</option>
    <option value="RU">Russian Federation</option>
    <option value="RW">Rwanda</option>
    <option value="BL">Saint Barthelemy</option>
    <option value="SH">Saint Helena</option>
    <option value="KN">Saint Kitts and Nevis</option>
    <option value="LC">Saint Lucia</option>
    <option value="MF">Saint Martin</option>
    <option value="PM">Saint Pierre and Miquelon</option>
    <option value="VC">Saint Vincent and the Grenadines</option>
    <option value="WS">Samoa</option>
    <option value="SM">San Marino</option>
    <option value="ST">Sao Tome and Principe</option>
    <option value="SA">Saudi Arabia</option>
    <option value="SN">Senegal</option>
    <option value="RS">Serbia</option>
    <option value="SC">Seychelles</option>
    <option value="SL">Sierra Leone</option>
    <option value="SG">Singapore</option>
    <option value="SK">Slovakia</option>
    <option value="SI">Slovenia</option>
    <option value="SB">Solomon Islands</option>
    <option value="SO">Somalia</option>
    <option value="ZA">South Africa</option>
    <option value="GS">South Georgia and the South Sandwich Islands</option>
    <option value="ES">Spain</option>
    <option value="LK">Sri Lanka</option>
    <option value="SD">Sudan</option>
    <option value="SR">Suriname</option>
    <option value="SJ">Svalbard and Jan Mayen</option>
    <option value="SZ">Swaziland</option>
    <option value="SE">Sweden</option>
    <option value="CH">Switzerland</option>
    <option value="SY">Syrian Arab Republic</option>
    <option value="TW">Taiwan</option>
    <option value="TJ">Tajikistan</option>
    <option value="TZ">Tanzania, United Republic of</option>
    <option value="TH">Thailand</option>
    <option value="TL">Timor-Leste</option>
    <option value="TG">Togo</option>
    <option value="TK">Tokelau</option>
    <option value="TO">Tonga</option>
    <option value="TT">Trinidad and Tobago</option>
    <option value="TN">Tunisia</option>
    <option value="TR">Turkey</option>
    <option value="TM">Turkmenistan</option>
    <option value="TC">Turks and Caicos Islands</option>
    <option value="TV">Tuvalu</option>
    <option value="UG">Uganda</option>
    <option value="UA">Ukraine</option>
    <option value="AE">United Arab Emirates</option>
    <option value="GB" selected>United Kingdom</option>
    <option value="US">United States</option>
    <option value="UM">United States Minor Outlying Islands</option>
    <option value="UY">Uruguay</option>
    <option value="UZ">Uzbekistan</option>
    <option value="VU">Vanuatu</option>
    <option value="VE">Venezuela</option>
    <option value="VN">Vietnam</option>
    <option value="VG">Virgin Islands, British</option>
    <option value="VI">Virgin Islands, U.S.</option>
    <option value="WF">Wallis and Futuna</option>
    <option value="EH">Western Sahara</option>
    <option value="YE">Yemen</option>
    <option value="ZM">Zambia</option>
    <option value="ZW">Zimbabwe</option>
</select>
</span>
        </label>
        <label class="input-set stacked input-set-state">
            <span class="label">State</span>
            <span class="input">
<select id="state-floating" name="XakFhrgoNgXPbP6AZtbKT1hVE3m8"  >
    <option value=""></option>
    <option value="AL">Alabama</option>
    <option value="AK">Alaska</option>
    <option value="AZ">Arizona</option>
    <option value="AR">Arkansas</option>
    <option value="CA">California</option>
    <option value="CO">Colorado</option>
    <option value="CT">Connecticut</option>
    <option value="DE">Delaware</option>
    <option value="DC">District of Columbia</option>
    <option value="FL">Florida</option>
    <option value="GA">Georgia</option>
    <option value="HI">Hawaii</option>
    <option value="ID">Idaho</option>
    <option value="IL">Illinois</option>
    <option value="IN">Indiana</option>
    <option value="IA">Iowa</option>
    <option value="KS">Kansas</option>
    <option value="KY">Kentucky</option>
    <option value="LA">Louisiana</option>
    <option value="ME">Maine</option>
    <option value="MD">Maryland</option>
    <option value="MA">Massachusetts</option>
    <option value="MI">Michigan</option>
    <option value="MN">Minnesota</option>
    <option value="MS">Mississippi</option>
    <option value="MO">Missouri</option>
    <option value="MT">Montana</option>
    <option value="NE">Nebraska</option>
    <option value="NV">Nevada</option>
    <option value="NH">New Hampshire</option>
    <option value="NJ">New Jersey</option>
    <option value="NM">New Mexico</option>
    <option value="NY">New York</option>
    <option value="NC">North Carolina</option>
    <option value="ND">North Dakota</option>
    <option value="OH">Ohio</option>
    <option value="OK">Oklahoma</option>
    <option value="OR">Oregon</option>
    <option value="PA">Pennsylvania</option>
    <option value="PR">Puerto Rico</option>
    <option value="RI">Rhode Island</option>
    <option value="SC">South Carolina</option>
    <option value="SD">South Dakota</option>
    <option value="TN">Tennessee</option>
    <option value="TX">Texas</option>
    <option value="UT">Utah</option>
    <option value="VT">Vermont</option>
    <option value="VA">Virginia</option>
    <option value="WA">Washington</option>
    <option value="WV">West Virginia</option>
    <option value="WI">Wisconsin</option>
    <option value="WY">Wyoming</option>
</select>
</span>
        </label>

        

<div class="gdpr-consent-requirement gdpr-consent-topics">
    <h4>
        Yes, also send me special offers about products &amp; services regarding:
        </h4>

        
        
 

<label class="input-set stacked inset">
    <span class="checkbox"> <input type="checkbox" name="XYl1zqRsaaqo5bEdEfbauXWzY5sg" value="artificial-intelligence"  data-consent-action data-consent-id=596517bdc14bed0737e41a50 ></span>
    <span class="checkbox-label">Artificial Intelligence</span>
    

</label>


        
        
 

<label class="input-set stacked inset">
    <span class="checkbox"> <input type="checkbox" name="XYl1zqRsaaqo5bEdEfbauXWzY5sg" value="cloud"  data-consent-action data-consent-id=596517bdc14bed0737e41a4f ></span>
    <span class="checkbox-label">Cloud</span>
    

</label>


        
        
 

<label class="input-set stacked inset">
    <span class="checkbox"> <input type="checkbox" name="XYl1zqRsaaqo5bEdEfbauXWzY5sg" value="network-security"  data-consent-action data-consent-id=596517bdc14bed0737e41a4e ></span>
    <span class="checkbox-label">Network Security</span>
    

</label>


        
        
 

<label class="input-set stacked inset">
    <span class="checkbox"> <input type="checkbox" name="XYl1zqRsaaqo5bEdEfbauXWzY5sg" value="hardware"  data-consent-action data-consent-id=596517bdc14bed0737e41a4d ></span>
    <span class="checkbox-label">Hardware</span>
    

</label>


        
        
 

<label class="input-set stacked inset">
    <span class="checkbox"> <input type="checkbox" name="XYl1zqRsaaqo5bEdEfbauXWzY5sg" value="software-development"  data-consent-action data-consent-id=596517bdc14bed0737e41a4c ></span>
    <span class="checkbox-label">Software Development</span>
    

</label>


        
</div>

<div class="gdpr-consent-requirement gdpr-contact-methods">
    <h4>
    You can contact me via:
    </h4>

    
 

<label class="input-set stacked inset input-set-consent-email minimum-explicit-required">
    <span class="checkbox"> <input type="checkbox" name="XYFF6oQoJeoESblJZdqpSqGJi58A" value="email"  data-consent-action data-consent-id=59aed8e456585fa9603b60ec ></span>
    <span class="checkbox-label">Email (required)</span>
    

</label>


    
 

<label class="input-set stacked inset prompt-phone">
    <span class="checkbox"> <input type="checkbox" name="XYFF6oQoJeoESblJZdqpSqGJi58A" value="phone"  data-consent-action data-consent-id=596517bdc14bed0737e41a51 ></span>
    <span class="checkbox-label">Phone</span>
    

</label>


    
 

<label class="input-set stacked inset prompt-phone">
    <span class="checkbox"> <input type="checkbox" name="XYFF6oQoJeoESblJZdqpSqGJi58A" value="sms"  data-consent-action data-consent-id=596517bdc14bed0737e41a53 ></span>
    <span class="checkbox-label">SMS</span>
    

</label>



    <label class="input-set inset input-phone hide ">
        <span class="label">Phone</span>
        <span class="input">
            <input type="tel" name="XakJ9oBANgXPbP6AZtbKT1hVE3m8" value="">
        </span>
        <span class="error-message"></span>
    </label>
</div>


        

      
        
    
    <div class="js-required fielderror">JavaScript is required for this form.</div>
    
    <div class="g-recaptcha"
          data-sitekey="6LeVgCEUAAAAACtawUTrPTBy0mTrGtjpPn_Xh-ZW"
          data-badge="inline"
          data-size="invisible">
    </div>


    </div>

    <p class="details">
    <span class="fielderror"></span>
    <input type="hidden" name="XZFxwuA0EfIE5ZkVehkOg523Q2os" value="sitewide" class="newsletter-optin-assume">
    
 

<label class="input-set input-set-agree-general allow-precheck">
    <span class="checkbox"> <input type="checkbox" name="XaFF6oQ0Nd4HVhaERifGdrF_pnlw" value="consent"  data-consent-action data-consent-id=596517bdc14bed0737e41a54 required></span>
    <span class="checkbox-label">I agree to receive correspondence from SourceForge.net.  I understand that I can withdraw my consent at anytime. Please refer to our <a href="http://slashdotmedia.com/terms-of-use">Terms of Use</a> and <a href="https://slashdotmedia.com/privacy-statement/">Privacy Policy</a> or <a href="/support">Contact Us</a> for more details.</span>
    

</label>


            
 

<label class="input-set input-set-agree-general-gdpr allow-precheck hide-initially">
    <span class="checkbox"> <input type="checkbox" name="XaFF6oQ0Nd4HVhaERifGdrF_pnlw" value="consent"  data-consent-action data-consent-id=596517bdc14bed0737e41a55 required></span>
    <span class="checkbox-label">I agree to receive correspondence from SourceForge.net via the means indicated above.  I understand that I can withdraw my consent at anytime. Please refer to our <a href="http://slashdotmedia.com/terms-of-use">Terms of Use</a> and <a href="https://slashdotmedia.com/privacy-statement/">Privacy Policy</a> or <a href="/support">Contact Us</a> for more details.</span>
    

</label>


    </p>

    <input type="submit" value="Subscribe" class="bt">

    <input type="hidden" name="source" value="floating">
    <input type="hidden" name="XfFF6uhAca4wSYFhJe5pFtI:WEpk" value="user">
    <input type="hidden" name="XbFl4uubuj8naN5xau8jZe1V3GAM" value="">


  <input type="hidden" name="_visit_cookie" value=""/>
<input type='hidden' name='timestamp' value='1505119436'/>
<input type='hidden' name='spinner' value='XbzIVz35oGfVNAzctHsU22_rkcfw'/>
<p class='H9effeddc4f6fef7a225a171e663950d2074984de'><label for='XaFp6oRsRKMXVhaERifGdrF_pnlw'>You seem to have CSS turned off.
             Please don't fill out this field.</label><input id='XaFp6oRsRKMXVhaERifGdrF_pnlw' name='XaVp6oRsRKW3LlQu6Km5BfvcL0VM' type=
             'text'/></p>
<p class='H9effeddc4f6fef7a225a171e663950d2074984de'><label for='XaFp6oRsRKMTVhaERifGdrF_pnlw'>You seem to have CSS turned off.
             Please don't fill out this field.</label><input id='XaFp6oRsRKMTVhaERifGdrF_pnlw' name='XaVp6oRsRKG3LlQu6Km5BfvcL0VM' type=
             'text'/></p>
</form>


    </div>
    <a id="btn-float-close">No, thanks</a>
</div>

        

<script src="//a.fsdn.com/con/js/sftheme/vendor/dropzone-4.3.0.min.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/vendor/dragster.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.metadata.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery-ui-1.11.1.custom.min.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.cookie.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.tablesorter.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.scrollTo-1.4.2.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/date.format.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/string.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/handlebars.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/lib/jquery.dotdotdot-1.8.3.min.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/jquery.notify.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/sftheme/shared.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/jquery.seltext.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/jquery.drawer.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/jquery.files.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/global.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/project.js?1504810730" type="text/javascript"></script>
<script src="//a.fsdn.com/con/js/madmen.js?1504810730" type="text/javascript"></script>



        
            <script type="text/javascript" src="//consent-st.truste.com/get?name=notice.js&amp;domain=slashdot.org&amp;c=teconsent&amp;text=true"></script>
            <noscript><p><img src="//analytics.slashdotmedia.com/sf.php?idsite=39" style="border:0;" alt="" /></p></noscript><script src="https://www.google.com/recaptcha/api.js?onload=recaptchaLoadCallback&render=explicit" async defer></script>
<script type="text/javascript" src="https://a.fsdn.com/con/js/files.js?1504810730"></script>


<script type="text/x-handlebars-template" id="file-drawer-template">

<form class="bp" action="{{files_url}}{{full_path}} method="put" id="file_properties_content">
    <table id="drawer_row">
        <col class="name-column">
        <col class="date-column">
        <col class="size-column">
        <col class="downloads-column">
        <col class="status-column">
        <tbody>
            <tr title="{{name}}">
                <td>
                {{#if authorized }}
                    <input type="text" class="title" name="name" value="{{name}}">
                {{else}}
                    <a href="{{file_title_url this}}" class="name">{{name}}</a>
                {{/if}}
                </td>
                <td class="files-date"></td>
                <td class="files-size"></td>
                <td class="files-downloads"></td>
                <td class="files-status status"></td>
            </tr>
            {{#if authorized}}
                <tr>
                    <td colspan="5" id="name_message" class="invalid hide"></td>
                </tr>
                {{#if d_type}}
                <tr>
                    <td colspan="5">
                        <input type="checkbox" name="stage" id="stage" value="1" style="vertical-align: baseline" {{checked this.staged}} {{stage_onclick this}} />
                        <label for="stage" title="Only release technicians will see this folder in the file browser." {{stage_onclick this}}>{{stage_message this staging_days}}</label>
                        <span title="Only release technicians will see this folder in the file browser." class="ico-help">?</span>
                        {{stage_date this}}
                    </td>
                </tr>
                {{/if}}
            {{/if}}
        </tbody>
    </table>

    <div id="file-details">
        {{#if f_type}}
        <div id="file-info">
            {{stage_date this}}

            <div class="label">
                <span>SHA1:</span>
            </div>
            <div class="value">
                <pre class="selectable">{{sha1}}</pre>
            </div>

            <div class="label">
                <span>MD5:</span>
            </div>
            <div class="value">
                <pre class="selectable">{{md5}}</pre>
            </div>

            {{#if authorized}}
            <div class="label">
                <span>Download URL:</span>
            </div>
            <div class="value">
                {{#if not_downloadable}}
                    <pre class="selectable" title="This file will be ready for download shortly.">This file will be ready for download shortly.</pre>
                {{else}}
                <pre class="selectable" title="{{download_url}}">{{download_url}}</pre>
                {{/if}}
            </div>

            <div class="label">
                <label for="download_label">Download Button:</label>
            </div>
            <div class="value">
                <span class="text"><input id="download_label" type="text" class="text" name="download_label" placeholder="Change the label used on the download button" value="{{download_label}}"></span>
            </div>

            <div class="label">
                <label for="exclude_reports">Exclude Stats:</label>
            </div>
            <div class="value">
                <span class="checkbox"><input type="checkbox" id="exclude_reports" name="exclude_reports" value="1" {{should_exclude_reports exclude_reports}}></span>
            </div>
            {{/if}}

            {{#if legacy_release_notes}}
            <div class="value no-label">
                <span><a href="{{legacy_release_notes}}">Release Notes</a></span>
            </div>
            {{/if}}
        </div>

        <div id="download-info">
            <div class="label" style="width: 137px">
                <span>Downloads (All-Time):</span>
            </div>
            <div class="value" style="width:100px">
                <span>{{download_display downloads}}</span>
            </div>

            {{#if authorized}}
            <div class="label">
                <span>Mirror Status:</span>
            </div>
            <div class="value">
                <span id="mirror_count">Loading ...</span>
            </div>
            {{/if}}

            <div class="default">
                {{#if show_platforms}}
                <div class="drawer-label">
                    <span>Default Download For:</span>
                </div>
                {{/if}}
                <div class="drawer-value">
                    <ul>
                    {{#each platforms}}
                        {{> platform}}
                    {{/each}}

                    {{#if authorized}}
                    <li><a href="#select_all" title="Select all">Select all</a></li>
                    {{/if}}
                    </ul>
                </div>
            </div>
        </div>
        {{/if}}

        {{#if authorized}}
        <hr />
        <p><input type="submit" value="Save"> <a href="#" id="cancel" class="btn link cancel">Cancel</a></p>
        {{/if}}
    </div>
</form>

</script>

<script type="text/x-handlebars-template" id="platform-partial">

    <li>
        <label>
            {{#if authorized}}
            <input type="checkbox" name="default" value="{{value}}" {{_checked}}>
            <span title="{{title}}" class="platform-icon {{value}}">{{title}}</span>
            {{/if}}

            {{#unless authorized}}
            {{#unless skip}}
            <span title="{{title}}" class="platform-icon {{value}}">{{title}}</span>
            {{/unless}}
            {{/unless}}
        </label>
    </li>

</script>

<script type="text/javascript">
    net.sf.files = {"fusioncatcher_v0.99.4b.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 207, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "ba5b5b2a6b59e4badba707148271b7519aabe9c4", "name": "fusioncatcher_v0.99.4b.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4b.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.4b.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "fb593bd2ad8713503f8dbf52dd44191d", "type": "f", "full_path": "fusioncatcher_v0.99.4b.zip"}, "fusioncatcher_v0.99.7a.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 61, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "0d7371c104a014197d57f54956f15e04208c45f9", "name": "fusioncatcher_v0.99.7a.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7a.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.7a.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "f50a50369d40ff1bd2f9bff29a35ef05", "type": "f", "full_path": "fusioncatcher_v0.99.7a.zip"}, "old": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 8, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "", "name": "old", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/old/download", "url": "/projects/fusioncatcher/files/old/", "downloadable": false, "authorized": null, "download_label": "", "md5": "", "type": "d", "full_path": "old"}, "extra": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 5, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "", "name": "extra", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/extra/download", "url": "/projects/fusioncatcher/files/extra/", "downloadable": false, "authorized": null, "download_label": "", "md5": "", "type": "d", "full_path": "extra"}, "examples": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 467, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "", "name": "examples", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/examples/download", "url": "/projects/fusioncatcher/files/examples/", "downloadable": false, "authorized": null, "download_label": "", "md5": "", "type": "d", "full_path": "examples"}, "fusioncatcher_v0.96.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 2, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "53288717c61c21809b9def5c42a7c587a277f60d", "name": "fusioncatcher_v0.96.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.96.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.96.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "78e5c95ecb16bb7a0e950a28c460dd5f", "type": "f", "full_path": "fusioncatcher_v0.96.zip"}, "fusioncatcher_v0.98.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 69, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "85bfa1fde7e8069d8f5a4b57c01c9adf1240bc17", "name": "fusioncatcher_v0.98.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.98.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.98.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "0fc50c93ef546463fd8db944f22bbaa9", "type": "f", "full_path": "fusioncatcher_v0.98.zip"}, "fusioncatcher_v0.99.7c.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 537, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "1ea6a1873cef12cd35ae0d7b6791934d313f9f7a", "name": "fusioncatcher_v0.99.7c.zip", "default": "linux", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7c.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.7c.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "1e7f71c6991c53c9b7ac9b185ba5f82e", "type": "f", "full_path": "fusioncatcher_v0.99.7c.zip"}, "fusioncatcher_v0.93.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 1, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "e62b0996185869394d44851c9eb2e413d7009e83", "name": "fusioncatcher_v0.93.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.93.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.93.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "a9cf977ce5ee8d45c706c8f47a0309f5", "type": "f", "full_path": "fusioncatcher_v0.93.zip"}, "fusioncatcher_v0.99.3a.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 215, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "3322778b0941d667d85efc7989e2a7a05ebbe662", "name": "fusioncatcher_v0.99.3a.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3a.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.3a.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "88136d5458bd693b074a7cf0ec1deb0e", "type": "f", "full_path": "fusioncatcher_v0.99.3a.zip"}, "fusioncatcher_v0.99.2.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 293, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "956a49dcfffeab41a80b9cde1e0b08c4eaaf7e88", "name": "fusioncatcher_v0.99.2.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.2.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.2.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "09171533b681bc2ed0179724f8d5291a", "type": "f", "full_path": "fusioncatcher_v0.99.2.zip"}, "fusioncatcher_v0.99.3c.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 173, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "6bf61c91ed6e6d0dade381a5a288c057bb57eca0", "name": "fusioncatcher_v0.99.3c.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3c.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.3c.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "3af9977c7fa881ceae4601efca9d521b", "type": "f", "full_path": "fusioncatcher_v0.99.3c.zip"}, "fusioncatcher_v0.94.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 3, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "6d28440e67b57fbb80e5217a2c61f80982c3e6d4", "name": "fusioncatcher_v0.94.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.94.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.94.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "24f07c886f23fca6dbef8f4a023e5578", "type": "f", "full_path": "fusioncatcher_v0.94.zip"}, "test": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 970, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "", "name": "test", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/test/download", "url": "/projects/fusioncatcher/files/test/", "downloadable": false, "authorized": null, "download_label": "", "md5": "", "type": "d", "full_path": "test"}, "fusioncatcher_v0.99.5a.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 571, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "6922b07a7097903529b356b0781c81f465fb5752", "name": "fusioncatcher_v0.99.5a.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.5a.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.5a.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "4d91c6afebb8b3979b8d2af6dceaca6f", "type": "f", "full_path": "fusioncatcher_v0.99.5a.zip"}, "fusioncatcher_v0.99.3b.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 390, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "d9b79c5e22216016266587f9f7a1da0e0d1ec4f1", "name": "fusioncatcher_v0.99.3b.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3b.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.3b.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "ef2b6a674d2c17f1c52da2655993f9da", "type": "f", "full_path": "fusioncatcher_v0.99.3b.zip"}, "fusioncatcher_v0.99.4d.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 500, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "981b733280f077a2c90870545efc42507d2523e8", "name": "fusioncatcher_v0.99.4d.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4d.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.4d.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "ce1b22a8067fd5ee2d30d72ce950be88", "type": "f", "full_path": "fusioncatcher_v0.99.4d.zip"}, "bootstrap.py": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 5219, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "b069f49287db90d51a7343520b742f92bfabdef4", "name": "bootstrap.py", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/bootstrap.py/download", "url": "/projects/fusioncatcher/files/bootstrap.py/", "downloadable": true, "authorized": null, "download_label": "", "md5": "8a1fdffe913f401db6bd4b087716d6d4", "type": "f", "full_path": "bootstrap.py"}, "fusioncatcher_v0.99.4c.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 386, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "ec5ea51257847636b1621fa02fe9d1646774b423", "name": "fusioncatcher_v0.99.4c.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4c.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.4c.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "130b4e2c467e0e3abeb515eddfff52a7", "type": "f", "full_path": "fusioncatcher_v0.99.4c.zip"}, "fusioncatcher_v0.99.6a.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 1730, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "89ae16811bfd21f9807119082cf983920d0793b2", "name": "fusioncatcher_v0.99.6a.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.6a.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "3b2306d35287626b3a4ca80d470e227f", "type": "f", "full_path": "fusioncatcher_v0.99.6a.zip"}, "data": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 9082, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "", "name": "data", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/data/download", "url": "/projects/fusioncatcher/files/data/", "downloadable": false, "authorized": null, "download_label": "", "md5": "", "type": "d", "full_path": "data"}, "README.md": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 47, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "1e93f502a0fdbce68ada6992a78818f1a4febea3", "name": "README.md", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/README.md/download", "url": "/projects/fusioncatcher/files/README.md/", "downloadable": true, "authorized": null, "download_label": "", "md5": "da1c8cc2b994f595fd606134d5087421", "type": "f", "full_path": "README.md"}, "fusioncatcher_v0.99.0.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 1, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "7e3e7986f370732e05a30b32c832dfa80a001581", "name": "fusioncatcher_v0.99.0.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.0.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.0.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "95abfb1e2272cc08811e2ef1bd1e41ff", "type": "f", "full_path": "fusioncatcher_v0.99.0.zip"}, "fusioncatcher_v0.99.3d.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 68, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "dd5c50adc660243424ad447887fe7baa6b869348", "name": "fusioncatcher_v0.99.3d.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3d.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.3d.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "1adc39bb0d552a2c9a8f88fc22014d12", "type": "f", "full_path": "fusioncatcher_v0.99.3d.zip"}, "fusioncatcher_v0.99.1.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 2, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "e8bdaacb7ea9c3bc71df3732294a1178648b4065", "name": "fusioncatcher_v0.99.1.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.1.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.1.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "89898a83a95b77f1a9e7dea570c73a98", "type": "f", "full_path": "fusioncatcher_v0.99.1.zip"}, "fusioncatcher_v0.97.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 4, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "5104b40283d7af10167fbf9b098185422c61c8fc", "name": "fusioncatcher_v0.97.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.97.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.97.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "a286d4ca33ca999b6609507097090da3", "type": "f", "full_path": "fusioncatcher_v0.97.zip"}, "fusioncatcher_v0.95.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 1, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "b8fbe404abbb3608039e61811170209f399f579c", "name": "fusioncatcher_v0.95.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.95.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.95.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "19d84aecd2333d48486a4d17f9cc072a", "type": "f", "full_path": "fusioncatcher_v0.95.zip"}, "fusioncatcher_v0.99.7b.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 683, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "e788c44c55f8099bc7cdc6a30c1d402351e6067d", "name": "fusioncatcher_v0.99.7b.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.7b.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.7b.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "1a35fe7a64cd1b6fc6e4b051f39fbc77", "type": "f", "full_path": "fusioncatcher_v0.99.7b.zip"}, "fusioncatcher_v0.99.3e.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 535, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "b7058be14cb41bede79c8e5ee5a98b3d8a064272", "name": "fusioncatcher_v0.99.3e.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.3e.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.3e.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "befac9fbfa263adb2f80e8a934f5e558", "type": "f", "full_path": "fusioncatcher_v0.99.3e.zip"}, "fusioncatcher_v0.99.4e.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 231, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "25ff5c0a0ba328d455f4ca905f1cd026a9eef047", "name": "fusioncatcher_v0.99.4e.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4e.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.4e.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "c1d31d79b0722965201a3d6a8f8c76ee", "type": "f", "full_path": "fusioncatcher_v0.99.4e.zip"}, "fusioncatcher_v0.99.4a.zip": {"staging_days": 3, "legacy_release_notes": null, "mirrors": 0, "downloads": 248, "exclude_reports": false, "explicitly_staged": false, "files_url": "/projects/fusioncatcher/files/", "link": "", "path": "", "stage": 0, "sha1": "d52948bcf4b4d098313b628ca932b019b2e8126a", "name": "fusioncatcher_v0.99.4a.zip", "default": "", "staged": false, "download_url": "https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v0.99.4a.zip/download", "url": "/projects/fusioncatcher/files/fusioncatcher_v0.99.4a.zip/", "downloadable": true, "authorized": null, "download_label": "", "md5": "9cb52358900bfc9b7f676d5bc3993199", "type": "f", "full_path": "fusioncatcher_v0.99.4a.zip"}};
    net.sf.staging_days = 3;
    $(function ($) {
        $('#files').files();
    });
    $(window).load(function() {
        $('.fs-sparklines').each(function(){
          $(this).attr('src', $(this).attr('data-src'));
        });

        var tsConfig = {
            widgets: ['zebra'],
            sortForce: [[4, 0]],
            headers: {
                1: { sorter: 'abbrdate' },
                2: { sorter: 'filesize' },
                3: { sorter: 'localnum' },
                4: { sorter: 'type'}
            },
            textExtraction: function (node) {
                return node.innerHTML;
            }
        };

        if ($('#files_list tbody tr').length) {
            $('#files_list').tablesorter(tsConfig);
            $('#files_list tr th#parent_folder').unbind();
            $('#files_list tr th.header').not('.typesort').not('#parent_folder').addClass('usersortable').append('<span class="ui-icon fs-toggle"><\/span>');
        }
    });
</script>

      </body>
</html>


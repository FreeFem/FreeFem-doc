const subnav = document.getElementById('subnav')
const contentdoc = document.getElementById('content-doc')
const contentexamples = document.getElementById('content-examples')
const examplewelcome = document.getElementById('exampleWelcome')
const examplecode = document.getElementById('exampleCode')
const examplelinks = document.getElementById('exampleLinks')
const tocdoc = document.getElementById('toc-doc')
const searchToggle = document.getElementById('checkbox-search-toggle')

let searchExamples = false

let viewMarkdown = true

let iscurrentfilemd = false

let mdsource = ""

toggleExamples = (mode) => {
  if (mode) {
    subnav.children[1].style.display = 'block'
    subnav.children[1].style.visibility = 'visible'
    subnav.children[0].style.display = 'none'
    subnav.children[0].style.visibility = 'hidden'

    contentexamples.style.display = 'block'
    contentexamples.style.visibility = 'visible'
    contentdoc.style.display = 'none'
    contentdoc.style.visibility = 'hidden'
    examplelinks.style.display = 'block'
    examplelinks.style.visibility = 'visible'
    tocdoc.style.display = 'none'
    tocdoc.style.visibility = 'hidden'
  } else {
    subnav.children[0].style.display = 'block'
    subnav.children[0].style.visibility = 'visible'
    subnav.children[1].style.display = 'none'
    subnav.children[1].style.visibility = 'hidden'

    contentdoc.style.display = 'block'
    contentdoc.style.visibility = 'visible'
    contentexamples.style.display = 'none'
    contentexamples.style.visibility = 'hidden'
    tocdoc.style.display = 'block'
    tocdoc.style.visibility = 'visible'
    examplelinks.style.display = 'none'
    examplelinks.style.visibility = 'hidden'
  }
}

loadExamplefromGitHub = (name, dir, editor) => {
  function loadedp(data) {
    iscurrentfilemd = false;
    var byClass = document.getElementsByClassName("strip-button");
    for (var i = 0; i < byClass.length; i++) {
      byClass[i].classList.add("disabled-strip-button");
    }
    var byClass = document.getElementsByClassName("button-strip");
    for (var i = 0; i < byClass.length; i++) {
      byClass[i].classList.add("button-strip-disabled");
    }
    examplecode.style.display = 'flex'
    examplewelcome.style.display = 'none'
    showExample(viewMarkdown)
    editor.setValue(data)
    highlightKeyword(editor)
  }

  function loadmd(data) {
    mdsource = data;
    iscurrentfilemd = true;
    var byClass = document.getElementsByClassName("strip-button");
    for (var i = 0; i < byClass.length; i++) {
      byClass[i].classList.remove("disabled-strip-button");
    }
    var byClass = document.getElementsByClassName("button-strip");
    for (var i = 0; i < byClass.length; i++) {
      byClass[i].classList.remove("button-strip-disabled");
    }
    //markdowntoggle.style.display = 'flex'
    examplecode.style.display = 'flex'
    examplewelcome.style.display = 'none'
    var edpcode = ""
    data.replace(/(```|~~~)freefem([\s\S]*?)\1/g, (match, p1,p2) => {
      edpcode += p2;
      return p2;
    });
    edpcode = edpcode.replace(/^\s/, '')
    editor.getWrapperElement().style.display = 'block'
    editor.getWrapperElement().style.visibility = 'hidden'
    editor.setValue(edpcode)
    highlightKeyword(editor)

    editor.getWrapperElement().style.display = 'none'
    editor.getWrapperElement().style.visibility = 'hidden'
    document.getElementById('mdout').innerHTML = md.render(data);

    MathJax.typeset();

    showExample(viewMarkdown);
  }

  const fpath = (dir == 'idp' ? '' : 'examples/') + dir + '/' + name
  const url = "https://raw.githubusercontent.com/FreeFem/FreeFem-sources/develop/" + fpath
  console.log('load ' + dir + '/' + name + ' from GitHub')
  document.getElementById('ExampleLinkToGitHub').innerHTML =
    "<a href='https://github.com/FreeFem/FreeFem-sources/blob/develop/" +
    fpath +
    "' target='_blank'>" +
    dir +
    '/' +
    name +
    '</a>'

  if (name.endsWith('.edp') || name.endsWith('.idp'))
    HTTPGet(url, loadedp)
  else
    HTTPGet(url, loadmd)
}

showExample = (viewMarkdown) => {
  if (viewMarkdown && iscurrentfilemd) {
    editor.getWrapperElement().style.display = 'none'
    editor.getWrapperElement().style.visibility = 'hidden'
    document.getElementById('mdout').style.display = 'block'
    document.getElementById('mdout').style.visibility = 'visible'
  }
  else {
    document.getElementById('mdout').style.display = 'none'
    document.getElementById('mdout').style.visibility = 'hidden'
    editor.getWrapperElement().style.display = 'block'
    editor.getWrapperElement().style.visibility = 'visible'
  }
}

FilterbyTag = (focus) => {
  filtered_examples.fill(true)
  const SetofTags = FilteredTags
  if (SetofTags.length != 0)
    all_examples.forEach(function callback(example, i) {
      const locTags = new Set(example.tags.split(' '))
      var take = true
      var n = 0

      while (n < SetofTags.length) {
        if (!locTags.has(SetofTags[n])) {
          take = false
          break
        }
        n++
      }

      if (!take) {
        filtered_examples[i] = false
      }
    })
}

updateList = (editor) => {
  var examplesList = document.getElementById('linksContainer')
  examplesList.innerHTML = ''
  var currentRep = ''
  all_examples.forEach(function callback(example, i) {
    if (filtered_examples[i] && filtered_examples_keyword[i]) {
      if (currentRep != example.dir) {
        examplesList.innerHTML += '<h3>' + example.dir + '/</h3>'
        currentRep = example.dir
      }

      const exi = document.createElement('a')
      exi.innerHTML =
        '<a href="#" title="example.desc" id="' +
        example.name +
        '" onclick="loadExamplefromGitHub(\'' +
        example.name +
        "','" +
        example.dir +
        '\', editor)"> ' +
        example.name +
        '</a><br>'
      examplesList.appendChild(exi)
    }
  })
}

highlightKeyword = (editor) => {
  let keyword = document.getElementById('searchInput').value
  if (keyword != '') {
    editor.doc.getAllMarks().forEach((marker) => marker.clear())
    var cursor = editor.getSearchCursor(keyword)
    while (cursor.findNext()) {
      editor.markText(cursor.from(), cursor.to(), {
        className: 'highlight-keyword'
      })
    }
  }
}

parsetree = (data) => {
  var countall = 0
  for (var i = 0; i < data.length; i++) {
    var count = 0
    if ('children' in data[i]) count = parsetree(data[i].children)
    else count = TagCount.get(data[i].id)
    data[i].text += ' (' + count + ')'
    countall += count
  }
  return countall
}

let treedata = [
  {
    id: '0',
    text: 'Equations',
    children: [
      {
        id: '0-0',
        text: 'Fluid mechanics',
        children: [
          { id: 'Stokes', text: 'Stokes' },
          { id: 'Navier-Stokes', text: 'Navier-Stokes' }
        ]
      },
      {
        id: '0-1',
        text: 'Solid mechanics',
        children: [{ id: 'elasticity', text: 'Elasticity' }]
      },
      {
        id: '0-2',
        text: 'Wave propagation',
        children: [
          {
            id: '0-2-0',
            text: 'Frequency domain',
            children: [
              { id: 'Helmholtz', text: 'Helmholtz equation' },
              { id: 'Maxwell', text: "Maxwell's equations" }
            ]
          }
          /*
    { "id": "0-2-1", "text": "Time domain"
    }
    */
        ]
      }
      /*
  { "id": "diffusion", "text": "Diffusion"
  },
  */
    ]
  },
  {
    id: '1',
    text: 'Spatial dimension',
    children: [
      { id: '2D', text: '2D' },
      { id: '3D', text: '3D' }
    ]
  },
  {
    id: '2',
    text: 'Mesh adaptation',
    children: [
      { id: 'adaptmesh', text: 'adaptmesh' },
      { id: 'tetg', text: 'tetg' },
      { id: 'mmg', text: 'mmg' }
    ]
  }
]

document.addEventListener(
  'keydown',
  (event) => {
    if (event.altKey && event.code == "AltRight") {
      searchExamples = !searchExamples
      toggleExamples(searchExamples)
      if (document.getElementById('search-overlay').style.display != 'none')
        search_refresh()
      if (searchExamples) searchToggle.checked = true
      else searchToggle.checked = false
    }
  },
  false
)

searchToggle.addEventListener('change', function () {
  if (this.checked) searchExamples = true
  else searchExamples = false
  toggleExamples(searchExamples)
  if (document.getElementById('search-overlay').style.display != 'none')
    search_refresh()
})

var editor = CodeMirror.fromTextArea(document.getElementById('code'), {
  lineNumbers: true,
  mode: 'text/x-ff++src'
})

searchExamples = false
toggleExamples(searchExamples)

var FilteredTags = []
var all_examples = []
var filtered_examples = []
var filtered_examples_keyword = []
const TagCount = new Map()
fetch('/_static/json/all_examples.json')
  .then((response) => response.json())
  .then((json) => {
    all_examples = json
    all_examples.forEach(function callback(example) {
      example.tags.split(' ').forEach(function callback(value) {
        TagCount.set(value, (TagCount.get(value) ?? 0) + 1)
      })
    })
    parsetree(treedata)
    filtered_examples = new Array(all_examples.length).fill(true)
    filtered_examples_keyword = new Array(all_examples.length).fill(true)
    updateList(editor)
    let tree = new Tree('.tag-tree-container', {
      data: treedata,
      closeDepth: 4,
      loaded: function () {
        //this.values = ['0-0-0', '0-1-1'];
        //console.log(this.selectedNodes);
        //console.log(this.values);
        //this.disables = ['0-0-0', '0-0-1', '0-0-2']
      },
      onChange: function () {
        FilteredTags = this.values
        FilterbyTag()
        updateList(editor)
      }
    })
  })

var copybutton = document.createElement('button')
copybutton.className = 'copy-button'
copybutton.style.top = '0rem'
copybutton.innerHTML = '<span class="clipboard-message">.edp copied to clipboard</span><i class="far fa-clone"></i>'
copybutton.onclick = function (e) {
  const textarea = document.createElement('textarea')
  if (iscurrentfilemd && viewMarkdown) {
    copybutton.innerHTML = '<span class="clipboard-message">.md copied to clipboard</span><i class="far fa-clone"></i>'
    textarea.value = mdsource
  }
  else {
    copybutton.innerHTML = '<span class="clipboard-message">.edp copied to clipboard</span><i class="far fa-clone"></i>'
    textarea.value = editor.getValue()
  }
  textarea.setAttribute('readonly', '')
  textarea.style.position = 'absolute'
  textarea.style.left = '-9999px'
  document.body.appendChild(textarea)
  textarea.select()
  document.execCommand('copy')
  document.body.removeChild(textarea)

  this.children[0].classList.toggle('clipboard-message--active')
  setTimeout(() => {
    this.children[0].classList.remove('clipboard-message--active')
  }, 2000)
}
const examplecopybutton = document.getElementById('ExampleCopyButton')
examplecopybutton.appendChild(copybutton)

var bs1 = new ButtonStrip({
  id: 'buttonStrip1'
});
bs1.addButton('.md', true, 'click', function(){
  viewMarkdown = true
  showExample(viewMarkdown)
});
bs1.addButton('.edp', false, 'click', function(){
  viewMarkdown = false
  showExample(viewMarkdown)
});
bs1.append('#markdowntoggle');
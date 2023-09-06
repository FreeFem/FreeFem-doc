const searchInput = document.getElementById('searchInput')
const searchResults = document.getElementById('searchResults')
let rootPath = ''
const LUNR_LIMIT = 10
const index = lunr.Index.load(LUNR_DATA)

const initSearch = (path) => {
  rootPath = path
}

const search = (event) => {
  searchClean()

  if (!event.target) return

  if (!event.target.value) return

  const text = event.target.value.toLowerCase()

  console.log(text)

  if (text === '') return
  else searchLunr(text)
}

const searchClean = () => {
  searchResults.innerHTML = ''
}

const searchLunr = (text) => {
  var results = index.query(function (q) {
    // look for an exact match and apply a large positive boost
    text.split(lunr.tokenizer.separator).forEach(function (term) {
      q.term(term, { usePipeline: true, boost: 10000 })

      // look for terms that match the beginning of this queryTerm and apply a medium boost
      if (text.length > 2) q.term(term, {
              wildcard: lunr.Query.wildcard.LEADING | lunr.Query.wildcard.TRAILING,
              usePipeline: false,
              boost: 100 }
            );

      // look for terms that match with an edit distance of 2 and apply a small boost
      if (text.length > 4) q.term(term, { usePipeline: false, editDistance: 2, boost: 1 })
    })
  })
  results = results.slice(0,20)

  // Search in page data
  const results_i = []
  results.map((result) => {
    const ref = result.ref
    const index_i = lunr.Index.load(LUNR_PAGEDATA[ref])
    results_i[ref] = index_i.query(function (q) {
      text.split(lunr.tokenizer.separator).forEach(function (term) {
        q.term(term, { usePipeline: true, boost: 10000 })

        // look for terms that match the beginning of this queryTerm and apply a medium boost
        if (text.length > 2) q.term(term, {
                wildcard: lunr.Query.wildcard.LEADING | lunr.Query.wildcard.TRAILING,
                usePipeline: false,
                boost: 100 }
              );

        // look for terms that match with an edit distance of 2 and apply a small boost
        if (text.length > 4) q.term(term, { usePipeline: false, editDistance: 2, boost: 1 })
      })
    })
  })

  const resultsHTML = parseLunrResults(results, results_i, text)
  if (resultsHTML.length === 0) searchResults.innerHTML = '<p>No results</p>'
  else {
    for (let result of resultsHTML) {
      const div = document.createElement('div')
      div.className = 'search-result'

      const title = document.createElement('a')
      title.className = 'search-result-main-title'
      title.innerHTML = '<i class="fas fa-file-alt"></i>' + result.mainTitle
      title.href = result.link
      div.appendChild(title)

      for (let i = 0; i < Math.min(30,result.titles.length); i++) {
        const subdiv = document.createElement('div')
        subdiv.className = 'search-result-sub'

        const subTitle = document.createElement('a')
        subTitle.className = 'search-result-title'
        subTitle.href = result.links[i]

        const titleDiv = document.createElement('div')

        const titleText = document.createElement('p')
        titleText.className = 'search-result-title-text'
        titleText.innerHTML = result.titles[i]

        const preview = document.createElement('p')
        preview.className = 'search-result-preview'
        preview.innerHTML = result.previews[i]

        subTitle.appendChild(titleDiv)

        titleDiv.appendChild(titleText)
        titleDiv.appendChild(preview)

        subdiv.appendChild(subTitle)

        div.appendChild(subdiv)
      }

      searchResults.appendChild(div)
    }
    searchResults.style.display = 'block'
  }
}

function parseLunrResults(results, results_i, text) {
  const html = []
  results.forEach((result) => {
    const ref = result.ref
    const item = PREVIEW_DATA[ref]
    const mainTitle = item.title
    const link = rootPath + item.link

    const titles = []
    const previews = []
    const links = []

    results_i[ref].forEach((result_i) => {
      const ref_i = result_i.ref
      const item_i = PREVIEW_PAGEDATA[ref][ref_i]
      const title_i = item_i.title
      const link_i = item_i.link
      var preview_i = item_i.preview
      var newpreview = ""

      Object.keys(result_i.matchData.metadata).forEach(function (term, index) {
        const lowerCaseText = term.toLowerCase()

        if (index == 0) {
          const lowerCasePreview = preview_i.toLowerCase()
          const textIndex = lowerCasePreview.indexOf(lowerCaseText)
          preview_i = preview_i.slice(
            Math.max(0, textIndex - 124),
            Math.min(textIndex + 124, preview_i.length)
          )
        }

        const lowerCasePreview2 = preview_i.toLowerCase()
        const textIndex2 = lowerCasePreview2.indexOf(lowerCaseText)
        if (textIndex2 != -1) {
          newpreview = newpreview +
            preview_i.slice(0, textIndex2) +
            '<b>' +
            preview_i.slice(textIndex2, textIndex2 + term.length) +
            '</b>'
          preview_i = preview_i.slice(textIndex2 + term.length)
        }
      })

      titles.push(title_i)
      previews.push(newpreview+preview_i)
      links.push(rootPath + link_i)
    })

    html.push({
      mainTitle,
      link,
      titles,
      previews,
      links
    })
  })

  return html
}

searchInput.addEventListener('input', function (event) {
  search(event)
})

searchInput.addEventListener('focus', function (event) {
  document.getElementById('search-overlay').style.display = 'block'

  if (!searchResults.children.length) return

  event.stopPropagation()
  searchResults.style.display = 'block'
})

function removeOverlay() {
  document.getElementById('search-overlay').style.display = 'none'
}

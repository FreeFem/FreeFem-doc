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

  const text = event.target.value

  if (text === '') return
  else searchLunr(text)
}

const searchClean = () => {
  searchResults.innerHTML = ''
}

const searchLunr = (text) => {
  const results = index.search(text)

  // Search in page data
  const results_i = []
  results.map((result) => {
    const ref = result.ref
    const index_i = lunr.Index.load(LUNR_PAGEDATA[ref])
    results_i[ref] = index_i.search(text)
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

      for (let i = 0; i < result.titles.length; i++) {
        const subdiv = document.createElement('div')
        subdiv.className = 'search-result-sub'

        const subtitle = document.createElement('a')
        subtitle.className = 'search-result-title'
        subtitle.href = result.links[i]

        const titleDiv = document.createElement('div')

        const titleText = document.createElement('p')
        titleText.className = 'search-result-title-text'
        titleText.innerHTML = result.titles[i]

        const preview = document.createElement('p')
        preview.className = 'search-result-preview'
        preview.innerHTML = result.previews[i]

        title.appendChild(titleDiv)

        titleDiv.appendChild(titleText)
        titleDiv.appendChild(preview)

        subdiv.appendChild(title)

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
      const preview_i = item_i.preview

      const lowerCasePreview = preview_i.toLowerCase()
      const lowerCaseText = text.toLowerCase()
      const textIndex = lowerCasePreview.indexOf(lowerCaseText)

      const splittedPreview = preview_i.slice(
        Math.max(0, textIndex - 124),
        Math.min(textIndex + 124, preview_i.length)
      )

      const lowerCasePreview2 = splittedPreview.toLowerCase()
      const textIndex2 = lowerCasePreview2.indexOf(lowerCaseText)
      const finalPreview =
        splittedPreview.slice(0, textIndex2) +
        '<b>' +
        splittedPreview.slice(textIndex2, textIndex2 + text.length) +
        '</b>' +
        splittedPreview.slice(textIndex2 + text.length)

      titles.push(title_i)
      previews.push(finalPreview)
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
